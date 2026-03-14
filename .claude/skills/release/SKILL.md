---
name: release
description: Register a new Julia package version with JuliaRegistrator including release notes. TagBot then creates the tag and GitHub release automatically.
disable-model-invocation: true
allowed-tools: Read, Grep, Bash(git *), Bash(gh *)
argument-hint: [version]
---

# Release a new version of the Julia package

Register a new version with JuliaRegistrator and auto-generated release notes. TagBot will create the git tag and GitHub release after the registry PR merges.

**Important:** Do NOT create the git tag manually. TagBot must create it so that it also creates the GitHub release. If the tag already exists, TagBot will skip the release.

## Instructions

1. **Determine the version to release:**
   - If `$ARGUMENTS` is provided, use it as the version (e.g., `v0.3.0`). Ensure it has a `v` prefix.
   - If no argument is given, read the `version` field from `Project.toml` and use that.

2. **Pre-flight checks:**
   - Run `git status` to ensure the working tree is clean. Abort if there are uncommitted changes.
   - Run `git tag` to verify the version tag does not already exist. Abort if it does — TagBot may have already released it, or a manual tag will block TagBot from creating the release.
   - Verify that the `version` field in `Project.toml` matches the version being released (without the `v` prefix). Warn the user if they don't match and ask how to proceed.

3. **Generate release notes** by running `git log <previous-tag>..HEAD --oneline` and categorising commits:
   - **Breaking changes** — commits with `!:` or `refactor!:`
   - **Bug fixes** — commits starting with `fix:`
   - **Improvements** — everything else worth noting (`feat:`, `docs:`, `ci:`, `build:`, `refactor:`, `perf:`)
   - Skip `style:` and `chore:` commits unless they are significant.
   - Write the notes in concise bullet points referencing function names where relevant.

4. **Show the user** the version and the draft release notes. Ask for confirmation before proceeding.

5. **Invoke JuliaRegistrator** by commenting on HEAD:
   ```
   gh api repos/{owner}/{repo}/commits/{sha}/comments -f body='@JuliaRegistrator register

   Release notes:

   <release notes here>'
   ```

6. **Report the result:** Confirm the registrator was invoked and link to the commit comment. Remind the user that:
   - The registry PR typically auto-merges after 3 days (shorter for patch/minor bumps of existing packages).
   - TagBot will automatically create the git tag and GitHub release (with the release notes) after the registry PR merges.
