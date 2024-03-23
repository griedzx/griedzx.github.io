@echo off
cd E:\my_website\
git add .
git commit -m "Daily update $(date)"
git push origin source
