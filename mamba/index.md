# Mamba

mamba安装

推荐安装方式

```shell
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
bash Miniforge3-$(uname)-$(uname -m).sh
```

1. `curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"`：这条命令使用 curl 工具从 GitHub 上下载 Miniforge 的最新版本。`-L` 参数让 curl 能够处理重定向，`-O` 参数让 curl 使用 URL 的文件名保存文件。`$(uname)` 和 `$(uname -m)` 是 shell 命令，它们分别返回操作系统的名称和机器硬件名称（例如，Linux x86_64），这样可以下载适合当前系统的 Miniforge 版本。
2. `bash Miniforge3-$(uname)-$(uname -m).sh`：这条命令使用 bash shell 执行下载的 Miniforge 安装脚本。

