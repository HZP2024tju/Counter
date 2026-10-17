# !/bin/bash
cd ..
rm -rf .git .vscode
git init
git add .
git commit -m "2026.10.15"
git remote add nerd git@github.com:HZP2024tju/Counter.git
git push --set-upstream nerd master -f