# Fanny Discord Bot

A Discord bot integrated with Nextflow to automate bioinformatics pipelines on HPC clusters.

## Installation

1.  **Clone the repository:**
    ```bash
    git clone <your-repo-url>
    cd fanny-bot
    ```

2.  **Sync Dependencies:**
    This project uses Go modules. Initialize the environment:
    ```bash
    go mod tidy
    ```

3.  **Configure:**
    Open `fanny_bot.go` and replace `YOUR_BOT_TOKEN_HERE` with your actual bot token from the Discord Developer Portal.

4.  **Compilation:**
    Build the binary:
    ```bash
    go build
    ```

## Usage

Set the following environment variables:

* DISCORD_BOT_TOKEN - the token for the bot
* DISCORD_CHANNEL - the channel the bot will run in
* DLE - the path to the folder containing the conda environments
* DLP - the path to the folder containing the nextflow pipelines

* **Launch the bot:**
    ```bash
    ./fanny-bot
    ```

* **Commands:**
    * `!fanny help` - Display command documentation.
    * `!fanny patho <run_id> [include-mito]` - Execute the pathogen pipeline for a given run.

## Requirements

* `nextflow` must be in the system PATH of the user running the bot.
* Go 1.18+
* Permissions: Ensure the bot has **Attach Files** and **Add Reactions** permissions in the channel.