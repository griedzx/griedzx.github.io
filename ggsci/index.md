# R绘图颜色配置ggsci

`ggsci` 是一个 R 包，提供了一些科学期刊、科学软件和科学网站的颜色主题。这些颜色主题可以直接用于 `ggplot2` 图形。

以下是如何使用 `ggsci` 的一些基本示例：

1. 加载 `ggsci` 包：

```r
library(ggsci)
```

2. 使用 `ggsci` 提供的颜色主题：

```r
library(ggplot2)

# 使用 NPG 颜色主题
ggplot(iris, aes(x = Sepal.Length, y = Sepal.Width, color = Species)) +
  geom_point() +
  scale_color_npg()
```

在这个例子中，`scale_color_npg` 函数被用于设置 NPG（Nature Publishing Group）的颜色主题。

3. 使用 `ggsci` 提供的调色板：

```r
# 显示 NPG 调色板的颜色
pal_npg()(4)
```

在这个例子中，`pal_npg` 函数被用于创建一个 NPG 调色板的函数，然后这个函数被调用来生成 4 种颜色。

`ggsci` 包提供了多种颜色主题和调色板，包括 `npg`、`aaas`、`lancet`、`jco`、`ucscgb`、`d3` 和 `gdocs` 等。你可以根据你的需要选择合适的颜色主题和调色板。

对不起，我误解了你的问题。你是对的，`ggsci` 包提供了离散和连续两种类型的颜色映射。

对于离散型变量，你可以使用 `scale_color_*` 或 `scale_fill_*` 函数，其中 `*` 是你选择的颜色主题。例如，如果你选择 `npg` 颜色主题，你可以使用 `scale_color_npg` 或 `scale_fill_npg` 函数。

对于连续型变量，你可以使用 `scale_color_continuous_*` 或 `scale_fill_continuous_*` 函数，其中 `*` 是你选择的颜色主题。例如，如果你选择 `npg` 颜色主题，你可以使用 `scale_color_continuous_npg` 或 `scale_fill_continuous_npg` 函数。

以下是一些使用 `ggsci` 包的例子：

```r
library(ggplot2)
library(ggsci)

# 使用离散型 NPG 颜色主题
ggplot(iris, aes(x = Sepal.Length, y = Sepal.Width, color = Species)) +
  geom_point() +
  scale_color_npg()

# 使用连续型 NPG 颜色主题
ggplot(faithful, aes(x = eruptions, y = waiting, color = waiting)) +
  geom_point() +
  scale_color_continuous_npg()
```

在这些例子中，`scale_color_npg` 和 `scale_color_continuous_npg` 函数被用于设置颜色的映射。你可以根据你的需要选择合适的颜色主题和映射函数。

[ggsci | ggplot2的颜色标度拓展包（1）：科研配色风格 - 知乎 (zhihu.com)](https://zhuanlan.zhihu.com/p/432716313)

