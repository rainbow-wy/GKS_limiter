# gks2d-str
这是一个气体动理学格式（Gas-Kinetic Scheme, GKS）的求解器，支持二维结构化网格模拟，目的是方便有需求的科研工作者入门之用。
GKS是香港科技大学华人学者徐昆教授原创的针对可压缩流动的格式（推荐一篇集大成的文献是 K. Xu, 2001, JCP, A Gas-Kinetic BGK Scheme for the Navier–Stokes Equations and Its Connection with Artificial Dissipation and Godunov Method）
如用作学术用途，请引用任意下列文献：
1) X. JI, F. ZHAO, W. SHYY, & K. XU (2018). A family of high-order gas-kinetic schemes and its comparison with Riemann solver based high-order methods. Journal of Computational Physics, 356, 150-173. 
2) X. JI, & K. XU (2020). Performance Enhancement for High-order Gas-kinetic Scheme Based on WENO-adaptive-order Reconstruction. Communication in Computational Physics, 28, 2, 539-590
## 功能
1. 一维网格 2阶到5阶 GKS
2. 二维均匀结构化单块网格 2阶到5阶 GKS
3. 二维非均匀结构化单块网格 2阶到5阶 GKS
4. 有限体积框架；WENO重构；多步多阶，两步四阶
## 安装
### windows 和 visual studio 2019
1. 创建新解决方案。 选择 “c++ windows桌面向导”，勾选 “将解决方案和项目放在同一目录”; 应用程序类型是 “控制台应用程序”; 选择 “空项目”
2. 把code文件夹（里面是所有源代码）和structured-mesh文件夹(是示例的结构化网格文件）复制到解决方案文件夹（和sln同目录）中。复制到其他文件夹或许也可以，关键点在于code文件夹和structured-mesh文件夹在一个目录下
3. 在visual studio中添加code中的 .h文件 到头文件；添加code中的.cpp文件到 源文件
4. 点击调试下拉菜单的 开始执行不调试（建议使用release（因为快）和x64（因为支持3.25GB以上的程序运行内存）编译）。编译成功，弹出控制台运行
### linux 和 g++
0. 目前linux下的安装教程可能是不完善的，因为我主要在windows下写程序
1. 复制code文件夹和和structured-mesh文件夹到同一目录下，确保安装了g++，支持c++11或者以上标准，安装了openmp
2. cd到code文件夹下目录，使用cmake生成makefile（cmake已经在code文件夹里面了），这一步是在code文件夹同级生成makefile
                                         cmake -B ../build
3. 使用cmake调用makefile生成可执行文件，这一步code文件夹同级生成build文件夹,可执行文件在build文件夹内
                                         cmake --build ../build
4. cd到build文件夹 执行可执行文件 main                             
                      
## 注                                        
程序 运行代码部分采用英文命名；注释部分采用中文+英文，以行文流畅易懂为准则，因此考虑到linux下不乱码，采用utf-8编码。
在windows系统VS下添加注释建议安装[forceUTF-8插件](https://blog.csdn.net/weixin_43734095/article/details/104481436) 为了兼顾vs下使用utf-8编码的注释不会扰乱运行代码，注释统一用//，而不用/**/              
20230720 bug 修正 vs2019下如果编译不通过可能是红山将文本文件的换行符改为\n(正则表达式)，使用 ctrl+f 搜索 \n （一定勾选使用正则表达式）替换为 \r\n 即可以符合windows的回车换行规则了。
# 致谢
感谢叶宇晨（北京大学）帮我测试了最初版本的代码  
感谢何志强博士告知我用git上传，尽管我现在仍一头雾水，只能用https 而不能用 ssh， 这是为什么？                              
特别鸣谢徐昆老师。该开源项目以及气体动理学讲习班（2022）是在徐昆老师的鼓励和支持下启动的。
                                                                                                    
