# Fix GitHub Tarball Symlink for Downstream pak

## Problem

The latest `prolfquapp` pkgdown workflow fails in `r-lib/actions/setup-r-dependencies@v2` while installing
`github::fgcz/prolfqua`. The failure happens before `prolfquapp` builds:

```text
Failed to uncompress prolfqua from ... prolfqua_1.6.1_02c6cf2.tar.gz-t
```

The GitHub tarball for `fgcz/prolfqua@02c6cf2` contains `AGENTS.md` as a symbolic link to `CLAUDE.md`.

## Fix

Replace the repository symlink `AGENTS.md -> CLAUDE.md` with a regular Markdown file. This keeps the same guidance text
but removes the symlink from GitHub-generated source archives consumed by `pak`.

## Validation

- Confirm `AGENTS.md` is a regular file in git.
- Confirm the package still parses/checks locally as far as practical.
- After pushing, rerun or wait for the downstream `prolfquapp` pkgdown workflow.
