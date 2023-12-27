# Azure_VPS

整Azure教育优惠，为期一年的超级自制VPS之旅 开整！！！

创建linux虚拟机


区域是服务器所在的数据中心位置，这个选项会影响连接服务器的延时。一般来说，选择周边的国家和地区速度会快一些。可以在[Azure Latency Test](https://link.zhihu.com/?target=https%3A//www.azurespeed.com/Azure/Latency)这个网站测试你的网络到不同数据中心的延时，也可以看到不同区域选项对应的真实地理位置，比如说East Asia对应的就是香港。如果你想通过服务器使用一些国外的服务，你也可能需要考虑一下服务器的位置。比如说要使用OpenAI的Chatgpt等服务，就不要选择East Asia。

![1703609642320](image/index/1703609642320.png)

映像就是服务器的操作系统。Linux系统的发行版众多，主要有Redhat和Debian两大分支。Redhat系下有RHEL，CentOS，Fedara等，Debian系下有Debian和Ubuntu等，除了这两系外还有许多其他的发行版。不同发行版的命令有一些差异，软件包也不同。映像可以根据你自己的需求选择。我看到网上相关教程比较多的主要是CentOS和Ubuntu这两个版本，这里选择Ubuntu Server 20.04作为示例

