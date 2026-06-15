/*
 * Copyright 2026 Christian Diener
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package main

import (
	"fmt"
	"log"
	"os"
	"os/exec"
	"os/signal"
	"path/filepath"
	"strconv"
	"strings"
	"sync"
	"syscall"

	"github.com/bwmarrin/discordgo"
)

var (
	BotToken      = os.Getenv("DISCORD_BOT_TOKEN")
	TargetChannel = os.Getenv("DISCORD_CHANNEL")
	Envs          = os.Getenv("DLE")
	Pipelines     = os.Getenv("DLP")
)

// pipelineLock ensures exclusive execution of the pipeline to prevent resource conflicts.
var pipelineLock sync.Mutex

// cleaner for input path
var replacer = strings.NewReplacer("`", "", "\"", "", "'", "")

func main() {
	dg, err := discordgo.New("Bot " + BotToken)
	if err != nil {
		log.Fatalf("Error creating session: %v", err)
	}

	dg.Identify.Intents = discordgo.IntentGuilds | discordgo.IntentsGuildMessages | discordgo.IntentsMessageContent
	dg.AddHandler(messageCreate)

	if err := dg.Open(); err != nil {
		log.Fatalf("Error opening connection: %v", err)
	}

	log.Println("Fanny Bot is running. Press CTRL-C to exit.")
	sc := make(chan os.Signal, 1)
	signal.Notify(sc, syscall.SIGINT, syscall.SIGTERM, os.Interrupt)
	<-sc
	dg.Close()
}

// messageCreate acts as the router for all incoming messages.
func messageCreate(s *discordgo.Session, m *discordgo.MessageCreate) {
	if m.Author.ID == s.State.User.ID {
		return
	}

	channel, _ := s.Channel(m.ChannelID)
	if channel == nil || channel.Name != TargetChannel {
		return
	}

	args := strings.Fields(m.Content)
	if len(args) == 0 || args[0] != "!fanny" {
		return
	}

	// Route to subcommands
	switch args[1] {
	case "help":
		sendHelp(s, m.ChannelID)
	case "talk":
		s.MessageReactionAdd(m.ChannelID, m.ID, "✅")
		s.ChannelMessageSend(m.ChannelID, fmt.Sprintf("👋 Hello %s, I am here and ready to help. 🧫", m.Author))
	case "patho":
		handlePatho(s, m, args[2:])
	default:
		s.ChannelMessageSend(m.ChannelID, "❌ Unknown subcommand. Use `!fanny help`.")
	}
}

// handlePatho manages the pathogen pipeline execution and result reporting.
func handlePatho(s *discordgo.Session, m *discordgo.MessageCreate, args []string) {
	if !pipelineLock.TryLock() {
		s.ChannelMessageSend(m.ChannelID, "⏳ **Pipeline is currently busy.** Please wait for the current run to finish.")
		return
	}
	defer pipelineLock.Unlock()

	if len(args) < 1 {
		s.ChannelMessageSend(m.ChannelID, "Usage: `!fanny patho <run_id> [include-mito]`")
		return
	}

	runArg := replacer.Replace(args[0])
	truncLen := 280
	if len(args) > 1 && args[1] == "include-mito" {
		truncLen = 205
	}

	folderDate := strings.ReplaceAll(strings.Split(runArg, "__")[0], "-", "")
	log.Printf(
		"Received a request for the Patho Pipeline for run %s with amplicon size %d. Will save results in %s.",
		runArg, truncLen, folderDate,
	)

	// Execute Nextflow Patho pipeline
	cmd := exec.Command(
		"nextflow", "run", "-resume",
		filepath.Join(Pipelines, "16S", "patho.nf"),
		"--run", runArg, "--read_length", strconv.Itoa(truncLen),
		"-with-conda", filepath.Join(Envs, "16S"),
		"-profile", "standard,collab",
	)

	s.ChannelMessageSend(
		m.ChannelID,
		fmt.Sprintf(
			"⏳ **The pipeline is now running** for run '%s'. \n"+
				"You can track the progress on https://cloud.seqera.io/orgs/dienerlab/workspaces/collab/watch.",
			runArg,
		),
	)

	if err := cmd.Run(); err != nil {
		s.MessageReactionAdd(m.ChannelID, m.ID, "❌")
		s.ChannelMessageSend(
			m.ChannelID,
			"⚠️ **Pipeline execution failed.**\n"+
				"Please check the logs and let Christian know in case it is not something obvious.",
		)
		return
	}
	s.MessageReactionAdd(m.ChannelID, m.ID, "✅")

	// Collect output files
	var files []*discordgo.File
	paths := []struct{ name, path string }{
		{"genus.png", filepath.Join(folderDate, "figures", "genus.png")},
		{"read_stats.csv", filepath.Join(folderDate, "tables", "read_stats.csv")},
	}
	for _, p := range paths {
		if f, err := os.Open(p.path); err == nil {
			files = append(files, &discordgo.File{Name: p.name, Reader: f})
			defer f.Close()
		} else {
			log.Printf("Could not read the file %s.", p)
		}
	}

	s.ChannelMessageSendComplex(m.ChannelID, &discordgo.MessageSend{
		Embeds: []*discordgo.MessageEmbed{{
			Title: "🚀 Pathogen Pipeline Finished",
			Description: fmt.Sprintf(
				`✅ **Patho Pipeline finished.**
I finished the pipeline and everything looks good.

You should be able to find the rest of the results in the folder '%s' at https://box.medunigraz.at/f/186637353 .
Attached are the logs, QC results and genus abundances.`,
				runArg,
			),
			Color: 0x00FF00,
			Fields: []*discordgo.MessageEmbedField{
				{Name: "Logs", Value: "```text\n" + getLogs(folderDate) + "\n```"},
			},
		}},
		Files: files,
	})
}

// getLogs aggregates and truncates log files for the result message.
func getLogs(folderDate string) string {
	files, _ := filepath.Glob(filepath.Join(folderDate, "logs", "*.log"))
	var sb strings.Builder
	for _, f := range files {
		if content, err := os.ReadFile(f); err == nil {
			sb.WriteString(fmt.Sprintf("--- %s ---\n%s\n", filepath.Base(f), string(content)))
		} else {
			log.Printf("Could not read the log %s.", f)
		}
	}
	res := sb.String()
	if len(res) > 1000 {
		return res[:1000] + "\n...[Truncated]"
	}
	return res
}

func sendHelp(s *discordgo.Session, cid string) {
	s.ChannelMessageSendEmbed(cid, &discordgo.MessageEmbed{
		Title: "🧬 Fanny Bot Help",
		Description: "I am Discord Bot that can run the our pipelines on an HPC.\n" +
			"I am named after Fanny Hesse and can be called with '!fanny' " +
			"followed by the subcommands listed below.",
		Fields: []*discordgo.MessageEmbedField{
			{
				Name:  "!fanny help",
				Value: "Show the available commands.",
			},
			{
				Name:  "!fanny talk",
				Value: "Check whether Fanny is there (the bot is connected).",
			},
			{
				Name: "!fanny patho <run_id> [include-mito]",
				Value: "Run the pathology 16S pipeline.\n" +
					"When the optional 'include-mito' part is provided it will choose a " +
					"shorter amplicon that will include mitochondria.",
			},
		},
		Color: 0x7289DA,
	})
}
