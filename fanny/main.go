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

const (
	CpuTotal   = 1416
	MemTotalGB = 9220.0
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

// Helper to determine emoji based on usage ratio
func getEmoji[T int | float64](used, total T) string {
	ratio := float64(used) / float64(total)
	switch {
	case ratio < 0.5:
		return "😇"
	case ratio < 0.8:
		return "🫡"
	case ratio <= 1.2:
		return "😰"
	default:
		return "😡"
	}
}

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
	case "cluster":
		if len(args) > 1 && args[2] == "status" {
			handleClusterStatus(s, m)
		} else {
			s.ChannelMessageSend(m.ChannelID, "❌ Unknown subcommand. Use `!fanny help`.")
		}
	case "patho":
		handlePatho(s, m, args[2:])
	default:
		s.ChannelMessageSend(m.ChannelID, "❌ Unknown subcommand. Use `!fanny help`.")
	}
}

// handleClusterStatus queries SLURM natively and safely.
func handleClusterStatus(s *discordgo.Session, m *discordgo.MessageCreate) {
	// Execute squeue without a shell. Includes both running and pending jobs.
	s.MessageReactionAdd(m.ChannelID, m.ID, "✅")
	cmd := exec.Command("squeue", "-h", "-t", "running,pending", "-p", "cpu", "-o", "%C %m")
	out, err := cmd.Output()
	if err != nil {
		log.Printf("Failed to query squeue: %v", err)
		s.ChannelMessageSend(m.ChannelID, "⚠️ Could not query cluster status.")
		return
	}

	cpuUsed := 0
	memUsedMB := 0.0

	lines := strings.Split(string(out), "\n")
	for _, line := range lines {
		fields := strings.Fields(line)
		if len(fields) < 2 {
			continue
		}

		// Parse CPUs
		c, _ := strconv.Atoi(fields[0])
		cpuUsed += c

		// Parse Memory (Handle suffixes if present)
		memStr := strings.TrimRight(fields[1], "M")
		m, _ := strconv.Atoi(memStr)
		memUsedMB += float64(m)
	}

	memUsedGB := memUsedMB / 1024.0

	_, err = s.ChannelMessageSendEmbed(m.ChannelID, &discordgo.MessageEmbed{
		Title:       "🖥️ MedBioNode status",
		Description: "Here is the status of the cluster (CPU partition).",
		Color:       0x7289DA,
		Fields: []*discordgo.MessageEmbedField{
			{
				Name:   "CPUs",
				Value:  fmt.Sprintf("%d/%d cores used %s", cpuUsed, CpuTotal, getEmoji(cpuUsed, CpuTotal)),
				Inline: true,
			},
			{
				Name:   "Memory",
				Value:  fmt.Sprintf("%.1f/%.1f GB used %s", memUsedGB, MemTotalGB, getEmoji(memUsedGB, MemTotalGB)),
				Inline: true,
			},
		},
	})

	if err != nil {
		log.Fatal("Failed to send the cluster status message.")
	}
}

// handlePatho manages the pathogen pipeline execution and result reporting.
func handlePatho(s *discordgo.Session, m *discordgo.MessageCreate, args []string) {
	if !pipelineLock.TryLock() {
		s.ChannelMessageSend(m.ChannelID, "⏳ **Patho Pipeline is currently busy.** Please wait for the current run to finish.")
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
			"⏳ **The pipeline is now running** for run `%s`. \n"+
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
		{"genus.png", filepath.Join(folderDate, "figures", "genus_top12.png")},
		{"read_stats.csv", filepath.Join(folderDate, "tables", "read_stats.csv")},
	}
	for _, p := range paths {
		if f, err := os.Open(p.path); err == nil {
			files = append(files, &discordgo.File{Name: p.name, Reader: f})
			defer f.Close()
		} else {
			log.Printf("Could not read the file %s.", p.path)
		}
	}

	_, err := s.ChannelMessageSendComplex(m.ChannelID, &discordgo.MessageSend{
		Embeds: []*discordgo.MessageEmbed{{
			Title: "🚀 Patho Pipeline Finished",
			Description: fmt.Sprintf(
				`✅ I finished the pipeline and everything **looks good**.

You should be able to find the rest of the results in the folder '%s' at https://box.medunigraz.at/f/186637353 .
Attached are the QC results and genus abundances.`,
				runArg,
			),
			Color: 0x00FF00,
		}},
		Files: files,
	})

	if err != nil {
		log.Fatal("Failed to send the final report message.")
	}
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
