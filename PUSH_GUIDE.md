# Push this folder to GitHub

Target repo: `https://github.com/newway515/late-life-depression-sex-gap-DIF`

## Step 1 — create an EMPTY repo on GitHub
Go to https://github.com/new
- Repository name: `late-life-depression-sex-gap-DIF`
- Public
- **Do NOT** check "Add a README", "Add .gitignore", or "Choose a license"
  (this folder already has them; adding them on the site causes a push conflict)
- Click **Create repository**

## Step 2 — push (run in this folder)
Open a terminal / PowerShell in `D:\ccproj\ccProject\ccProj01\github_repo`, then:

```bash
git init
git add .
git commit -m "Analysis code for late-life depression sex-gap DIF study"
git branch -M main
git remote add origin https://github.com/newway515/late-life-depression-sex-gap-DIF.git
git push -u origin main
```

When prompted to sign in, a browser window will open — log in as `newway515` and authorize.

### Alternative: GitHub CLI (creates the repo AND pushes in one go)
If you have `gh` installed and logged in (`gh auth login`), skip Step 1 and run:
```bash
git init && git add . && git commit -m "Analysis code for late-life depression sex-gap DIF study"
gh repo create newway515/late-life-depression-sex-gap-DIF --public --source=. --remote=origin --push
```

## Step 3 — verify
Open https://github.com/newway515/late-life-depression-sex-gap-DIF and confirm you see
`R/`, `python/`, `README.md`, `LICENSE`. The manuscript's Code-availability link now resolves.

---

## Common snags
- **"git is not recognized"** → install Git for Windows (https://git-scm.com/download/win), reopen the terminal.
- **"remote origin already exists"** → run `git remote remove origin` then re-add.
- **Push rejected / "updates were rejected"** → the GitHub repo wasn't empty (you added a README/license on the site). Either recreate it empty, or run `git pull --rebase origin main` then push again.
- **Wrong files / data leaked?** The included `.gitignore` blocks `*.db`, `*.csv`, `*.rds`, `data/`, `export/`, `output/`. Before committing you can preview what will be tracked with `git status`; nothing under those patterns should appear.
- **Repo name** must stay `late-life-depression-sex-gap-DIF` to match the link already written into the manuscript. If you change it, tell me and I'll update the docx.
