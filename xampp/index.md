# xampp集成软件包配置

## XAMPP VPS(Linux)下使用

### 下载

[Download xampp-linux-x64-8.2.12-0-installer.run (XAMPP) (sourceforge.net)](https://sourceforge.net/projects/xampp/files/XAMPP%20Linux/8.2.12/xampp-linux-x64-8.2.12-0-installer.run/download)

没有找到直接下载的地址，在本地下载安装脚本，scp上传至服务器/VPS上

### 修改端口号

服务器上通常会有一些其他的常见web应用占用了端口号，会和xampp下软件默认端口冲突

1. Another web server is already running

   说明apache端口号被占用

   进入/opt/lampp/etc目录下，找到文件httpd.conf，此文件中有两个地方需要修改

   a.找到Listen行修改，如：1211

   b.找到ServerName行修改，如：10.0.0.251:1211
2. Another web server with SSL is already running

   说明SSL端口号被占用

   进入/opt/lampp/etc/extra目录下，找到文件httpd-ssl.conf，此文件中有三个地方需要修改

   a.找到Listen行修改，如：1212

   b.找到ServerName行修改，如：10.0.0.251:1212

   c.找到节点VirtualHost行修改，如：_default_:1212
3. Another MySQL daemon is already running

   说明MySQL端口号被占用

   进入/opt/lampp/etc目录下，找到文件my.cnf，此文件中有两个地方需要修改

   a.分别找到port行修改，如：1213
4. Another FTP daemon is already running

   说明FTP端口号被占用

   进入/opt/lampp/etc目录下，找到文件proftpd.conf，此文件中有一个地方需要修改

   a.找到Port行修改，如：1214

修改完上面各自的配置文件端口号后，再修改xampp中的配置端口号

进入/opt/lampp目录下，找到文件xampp，此文件中有四个地方需要修改
找到testport，总共有四个，分别对应apache、SSL、MySql、FTP的端口号根据上面修改的端口号，对应修改此配置文件中的端口号（上面没有修改的端口，则对应的不需要改动）

### 其他报错

```text
(base) [root@dzx lampp]# ./xampp start
Starting XAMPP for Linux 8.2.12-0...
XAMPP: Starting Apache...fail.
httpd: Syntax error on line 522 of /opt/lampp/etc/httpd.conf: Syntax error on line 6 of /opt/lampp/etc/extra/httpd-xampp.conf: Cannot load modules/mod_perl.so into server: libnsl.so.1: cannot open shared object file: No such file or directory
XAMPP: Starting MySQL...ok.
XAMPP: Starting ProFTPD...already running.
```

XAMPP启动失败，主要是因为Apache无法加载mod_perl.so模块，因为它依赖的libnsl.so.1库文件找不到，在某些Linux发行版中，libnsl.so.1库可能不再默认安装。`sudo yum install libnsl`

```text
2024-03-29 6:32:00 0 [ERROR] Can't start server: Bind on TCP/IP port. Got error: 98: Address already in use 
2024-03-29 6:32:00 0 [ERROR] Do you already have another mysqld server running on port: 3306 ? 
2024-03-29 6:32:00 0 [ERROR] Aborting 2024-03-29 06:32:01 51844 mysqld_safe mysqld from pid file /opt/lampp/var/mysql/dzx.pid ended 
2024-03-29 06:32:13 52226 mysqld_safe Starting mysqld daemon with databases from /opt/lampp/var/mysql
```

检查xampp 安装目录中 /opt/lampp/var/mysql/*.err 一个错误日志，在底部会有具体的问题信息，我的问题是映射到了一个其他的docker应用使用的mysql

[MySQL is not working in XAMPP 7.3.0-0 &#34;/opt/lampp/bin/mysql.server: 260: kill: No such process&#34; - Ask Ubuntu](https://askubuntu.com/questions/1102355/mysql-is-not-working-in-xampp-7-3-0-0-opt-lampp-bin-mysql-server-260-kill-n)

上述问题解决后运行 `xampp`还是有报错，但是不影响xampp使用

```text
(base) [root@dzx lampp]# /opt/lampp/bin/mysql.server: line 262: log_success_msg: command not found
```

可以看到是一个日志输出，/opt/lampp/bin/mysql.server脚本中一个命令没有找到

[php - 如何修复 /opt/lampp/bin/mysql.server：第 261 行：log_success_msg：找不到命令？- 堆栈溢出 (stackoverflow.com)](https://stackoverflow.com/questions/72753116/how-to-fix-opt-lampp-bin-mysql-server-line-261-log-success-msg-command-not-f)

可以直接删除对应函数，也可以更新下载对应命令  `sudo yum install redhat-lsb-core`

### 通过ip/域名访问报错

网页报错

```plaintext
griedzx: Access forbidden!

New XAMPP security concept:

Access to the requested directory is only available from the local network.

This setting can be configured in the file "httpd-xampp.conf".
```

这个错误信息表明，正在尝试从非本地网络访问XAMPP，但是XAMPP的安全设置不允许这样做。默认情况下，XAMPP只允许从本地网络访问。根据提示可以修改 `httpd-xampp.conf`文件来改变这个设置。这个文件通常位于 `/opt/lampp/etc/extra/`目录下。

在 `httpd-xampp.conf`文件中，找到以下部分：

通过vim可以在底层命令行通过 `/<text>`快速查找相关内容

```apache
<LocationMatch "^/(?i:(?:xampp|security|licenses|phpmyadmin|webalizer|server-status|server-info))">
    Require local
    ErrorDocument 403 /error/XAMPP_FORBIDDEN.html.var
</LocationMatch>
```

将 `Require local`改为 `Require all granted`，如下所示：

```apache
<LocationMatch "^/(?i:(?:xampp|security|licenses|phpmyadmin|webalizer|server-status|server-info))">
    Require all granted
    ErrorDocument 403 /error/XAMPP_FORBIDDEN.html.var
</LocationMatch>
```

保存文件并重启XAMPP，就可以从非本地网络访问XAMPP了。

这样做可能会增加安全风险，因为任何人都可以访问这个XAMPP

注意， 当访问一个目录（例如 `http://xampp.daizixi.space:1080/2021317210202/`）而不是具体的文件时，Web服务器会查看该目录下是否存在这些默认文件 `（index.html`，`index.php`等）。如果存在，服务器就会自动返回这个文件的内容

类似的，正在使用的hugo框架也是相应的规则 [hugo主题框架中文经验贴](https://www.jianshu.com/p/0b9aecff290c)


## XAMPP windows下使用

[XAMPP download | SourceForge.net](https://sourceforge.net/projects/xampp/)

window下本地使用类似，且没什么其他应用也不存在端口冲突，并且网上window下教程很多，不像linux，让我折腾了将近一天

## 参考

[linux下用xampp安装php集成环境，并修改各自端口号 - 刚子2013 - 博客园 (cnblogs.com)](https://www.cnblogs.com/gangzi2013/p/6030157.html)

[MySQL is not working in XAMPP 7.3.0-0 &#34;/opt/lampp/bin/mysql.server: 260: kill: No such process&#34; - Ask Ubuntu](https://askubuntu.com/questions/1102355/mysql-is-not-working-in-xampp-7-3-0-0-opt-lampp-bin-mysql-server-260-kill-n)

[php - How to fix /opt/lampp/bin/mysql.server: line 261: log_success_msg: command not found? - Stack Overflow](https://stackoverflow.com/questions/72753116/how-to-fix-opt-lampp-bin-mysql-server-line-261-log-success-msg-command-not-f)

