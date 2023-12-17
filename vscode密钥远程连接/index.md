# Vscode密钥远程连接

通过vscode的remote-ssh插件可以在vscode上远程登陆服务器，登录和打开工作文件夹时都需要输入密码，可以通过配置ssh免密登录，使登录和使用更加丝滑

## 必要软件安装

`ssh-keygen` 和 `ssh-copy-id` 安装，window电脑powershell中这两个命令已内置

## ssh-keygen配置密钥对

配置密钥

```powershell
PS C:\Users\griedzx> ssh-keygen
Generating public/private rsa key pair.
Enter file in which to save the key (C:\Users\griedzx/.ssh/id_rsa):
Enter passphrase (empty for no passphrase):
Enter same passphrase again:
Your identification has been saved in C:\Users\griedzx/.ssh/id_rsa
Your public key has been saved in C:\Users\griedzx/.ssh/id_rsa.pub
The key fingerprint is:
SHA256:sYURFrXR+/Xffdy+HmOSu1yxPCH8WbXcoKNvSZF4RDQ griedzx@dzx_pc
The key's randomart image is:
+---[RSA 3072]----+
|        =+o+E    |
|       . o oo.   |
|        o oo o. .|
|         +. *...=|
|        S  .o= *+|
|           ...* B|
|          .. + %=|
|           .+ = @|
|           ..+o++|
+----[SHA256]-----+
```

注意 `Enter passphrase (empty for no passphrase):` 不输入密码，直接回车继续，否则之前对应输入服务器用户密码的地方变成了密钥对的密码

ssh-keygen常用参数

```plaintext
ssh-keygen
-t: 密钥类型, 可以选择 dsa | ecdsa | ed25519 | rsa;

-f: 密钥目录位置, 默认为当前用户home路径下的.ssh隐藏目录, 也就是~/.ssh/, 同时默认密钥文件名以id_rsa开头. 如果是root用户, 则在/root/.ssh/id_rsa, 若为其他用户, 则在/home/username/.ssh/id_rsa;

-C: 指定此密钥的备注信息, 需要配置多个免密登录时, 建议携带;

-N: 指定此密钥对的密码, 如果指定此参数, 则命令执行过程中就不会出现交互确认密码的信息了.
```

设定目录下会出现密钥对文件：id_rsa.pub是公钥 id_rsa是私钥


## 公钥传输至服务器

公钥放server(远程主机)上，私钥放本机上

### 代码传输

公钥配置至服务器上（Window -> Linux）

```powershell
$USER_AT_HOST = "yh@760755rf68.imdo.co" #id@ip
$PUBKEYPATH = "$HOME/.ssh/id_rsa.pub"

$pubKey = (Get-Content "$PUBKEYPATH"|Out-string); ssh -p 37763 "$USER_AT_HOST" "echo '${pubKey}' >> ~/.ssh/authorized_keys"
```



### 简单粗暴的方法

建议直接cat 本地的id_rsa.pub

```powershell
PS C:\Users\griedzx\.ssh> cat id_rsa.pub
```



然后复制内容到服务器的~/.ssh/authorized_keys中新增一行

#服务器端

```shell
(base) yh@localhost 21:21:38 ~/.ssh
$ vim authorized_keys
```



## 修改$home/.ssh/config

修改.ssh/config文件：加入IdentityFile的路径（也就是私钥在本机的所在位置）

```powershell
Host yh_xw
  HostName 760755rf68.imdo.co
  Port 37763
  User yh
  IdentityFile "C:\Users\griedzx\.ssh\id_rsa"
```



vscode登录server就不用输入密码了

