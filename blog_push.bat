@echo off
setlocal enabledelayedexpansion

:: 定义输出文件名
set "outputFile=history.txt"

:: 打印日期并重定向到文件
echo %date% >> %outputFile%

:: 计算文件中的日期行数
set "lineCount=0"
for /f %%A in ('type "%outputFile%"^|find /c /v ""') do set "lineCount=%%A"

:: 如果日期行数达到5，覆盖文件
if !lineCount! geq 5 (
    echo %date% > %outputFile%
)

cd E:\my_website\
git add .
git commit -m "Daily update %date%"
git push origin source >> %outputFile% 2>&1