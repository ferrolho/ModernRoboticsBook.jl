---
name: release
description: Tag a new Julia package version, push the tag, and register with JuliaRegistrator including release notes
disable-model-invocation: true
allowed-tools: Read, Grep, Bash(git *), Bash(gh *)
argument-hint: [version]
---

# Release a new version of the Julia package

Tag a new version, push it, and register with JuliaRegistrator along with auto-generated release notes.

## Instructions

1. **Determine the version to release:**
   - If `$ARGUMENTS` is provided, use it as the version (e.g., `v0.3.0`). Ensure it has a `v` prefix.
   - If no argument is given, read the `version` field from `Project.toml` and use that.

2. **Pre-flight checks:**
   - Run `git status` to ensure the working tree is clean. Abort if there are uncommitted changes.
   - Run `git tag` to verify the version tag does not already exist. Abort if it does.
   - Verify that the `version` field in `Project.toml` matches the tag being created (without the `v` prefix). Warn the user if they don't match and ask how to proceed.

3. **Generate release notes** by running `git log <previous-tag>..HEAD --oneline` and categorising commits:
   - **Breaking changes** — commits with `!:` or `refactor!:`
   - **Bug fixes** — commits starting with `fix:`
   - **Improvements** — everything else worth noting (`feat:`, `docs:`, `ci:`, `build:`, `refactor:`, `perf:`)
   - Skip `style:` and `chore:` commits unless they are significant.
   - Write the notes in concise bullet points referencing function names where relevant.

4. **Show the user** the tag version and the draft release notes. Ask for confirmation before proceeding.

5. **Create and push the tag:**
   ```
   git tag <version>
   git push origin <version>
   ```

6. **Invoke JuliaRegistrator** by commenting on the tagged commit:
   ```
   gh api repos/{owner}/{repo}/commits/{sha}/comments -f body='@JuliaRegistrator register

   Release notes:

   <release notes here>'
   ```

7. **Report the result:** Confirm the tag was pushed and link to the commit comment. Remind the user that the registry PR typically auto-merges after 3 days, and TagBot will create the GitHub release.
