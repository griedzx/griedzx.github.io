# 搭建博客

# 搭建博客记录

## 选择Hugo&Githubpages的原因

1. Hugo是一个用Go语言编写的静态网站生成器，它可以快速地创建和发布博客文章。基础运用可以简单，只需要安装Hugo，选择一个主题，编写Markdown格式的文章，然后使用Hugo命令生成静态网站文件，不需要担心数据库、服务器等问题
2. GitHub Pages是一个静态站点托管服务，直接将个人、组织或项目的页面托管于GitHub库或仓库中。使用GitHub可以创建一个 `<username>.github.io`的网站，不需要第一时间考虑域名以及备案等问题，后续有时间可以折腾--更改至自己的域名
3. GitHub workflow 可以自动更新Hugo 博客，这是一个很方便的功能。只需要将 Hugo 源码 push 到 GitHub 上，就可以触发 GitHub Actions 来生成和部署您的静态网站文件。您不需要在本地运行 Hugo 命令，也不需要手动上传文件到服务器
4. LoveIt主题个人比较喜欢，相对wordpress上大部分跟简洁，想要花哨也可以自定义魔改主题👴

## 基于Hugo&Githubpages的博客搭建流程

### 基本流程

1. 安装 Hugo。可以从[Hugo官网]下载适合您的操作系统的版本，建议下载entende版本（**由于很多主题的一些特性需要将 SCSS 转换为 CSS, 推荐使用 Hugo extended版本来获得更好的使用体验** ）
   可以使用 `hugo version` or `hugo -h`来检验是否安装成功
2. 创建一个新的 Hugo 站点。可以使用 `hugo new site path/to/site`命令在指定的路径下创建一个空的站点框架
3. 下载主题模板，我使用的是Git 模块方式添加主题到指定文件夹下 `git submodule add https://github.com/griedzx/LoveIt.git themes/LoveIt`
   我将LoveIt模板文件fork到自己的仓库，并使用 `git submodule`（日后自己个性化模板，可以一键同步更新）
4. 相应设置 config.toml，指定主题、网址，可以参考对应主题中的toml文件
5. `hugo new posts/<xxx>.md`创建markdown文件作文博文
6. `hugo serve -D` 启动 Hugo server 并使用 drafts 模式，在本地运行一个 web 服务器来预览网站效果，并且可以看到草稿文件的预览
7. 使用GitHub可以创建一个 `<username>.github.io`的仓库，将 `hugo`后产生的 `./public`（里面存放的是可部署的静态网站文件）push至对应repo
   7.1 GitHub workflow 可以自动更新您的 Hugo 博客，这是一个很方便的功能。只需要将您Hugo 源码 push 到 GitHub 上（可以是另一个repo，也可以是上述repo另一个分支），就可以触发 GitHub Actions 来自动生成和部署的静态网站文件

   如果想了解更多关于 GitHub workflow 的信息，您可以参考[Host on Github Pages | Hugo](https://gohugo.io/hosting-and-deployment/hosting-on-github)

### 设置的tips

#### 文件头扉页 front matter

```
title: "搭建博客"
date: 2023-07-28T21:57:35+08:00
draft: false
```

hugo一个新md文件一般只会出现前三行设置参数，但是Hugo 将内容分成草稿 Draft，将来发布 Future 和过期 Expired 等类型，可以在文件头扉页 front matter 中设置相应状态。

* future 设置 publishdate 值
* draft 设置 true 或者 false
* past 设置 expirydate 值
* 如 demo.md 文件头扉页 front matter 中设置：

```

---
title: Base Templates and Blocks
linktitle:
description: The base and block constructs ...
godocref: https://golang.org/pkg/text/template/#example_Template_block
date: 2017-02-01
publishdate: 2017-02-01
lastmod: 2017-02-01
categories: [templates,fundamentals]
keywords: [blocks,base]
menu:
  docs:
    parent: "templates"
    weight: 20
weight: 20
sections_weight: 20
draft: false
aliases: [/templates/blocks/,/templates/base-templates-and-blocks/]
toc: true
---
```



#### 多语言模式

[language]下设置多种参数，具体见[多语言模式 |雨 果 (gohugo.io)](https://gohugo.io/content-management/multilingual/)

对于每个新页面, 将语言代码附加到文件名中.

单个文件 `my-page.md` 需要分为三个文件:

* 英语: `my-page.en.md`
* 中文: `my-page.zh-cn.md`
* 法语: `my-page.fr.md`

#### emoji展示

win10输入法里面 徽标 + 句号（.）可以使用自带的表情符号，可以个性化设置 🐱‍🏍 好玩



以上完成初步hugo网站的搭建，还有很多设置可以学习。hugo官方的汉化版使用文档很简略，大部分网上找到的翻译文档大多直译看完后还是摸不着头脑，在简书上发现[Hugo 不完美教程 ](https://www.jianshu.com/p/deaa0e58315a)这一系列教程，使用文档和作者实战经验相结合，对我有很大帮助，我将借助此继续完善我的网站


参考

[Hugo 不完美教程 - I: Hugo Web Framework - 简书 (jianshu.com)](https://www.jianshu.com/p/deaa0e58315a)

[Hugo框架中文文档 标签分类 - Andbible](https://www.andbible.com/post/hugo-content-management-taxonomies/)

