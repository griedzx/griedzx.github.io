@echo off
setlocal enabledelayedexpansion

:: 定义输出文件名
set "outputFile=E:\my_website\static\blog.log"

:: 打印日期并重定向到文件
echo %date% >> %outputFile%

:: 计算文件中的日期行数
set "dateCount=0"
for /f "delims=" %%A in ('findstr /r "^[0-9][0-9]*-[0-9][0-9]*-[0-9][0-9]*" "%outputFile%"') do (
    set /a "dateCount+=1"
)

:: 如果日期行数达到5，覆盖文件
if !dateCount! geq 5 (
    echo %date% > %outputFile%
)

cd E:\my_website\ >> %outputFile% 2>&1
git add . >> %outputFile% 2>&1
git commit -m "Daily update %date%" >> %outputFile% 2>&1
git push origin source >> %outputFile% 