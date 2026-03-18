# HealShapes

### 本人环境

- IDE：VS2022

  - 语言标准：C++20

    - 新特性用的不是很多，就用了`std::format`、`std::string_view`，看情况是否降级

  - 字符集：UTF-8，记得在命令行中添加对应指令（因为师兄用的CMake + GCC，所以是UTF-8字符集）

  - 附加依赖项：

    ```
    geom.lib;topo.lib;nlopt.lib;TKernel.lib;TKMath.lib;TKBRep.lib;TKG2d.lib;TKG3d.lib;TKGeomBase.lib;TKDESTEP.lib;TKDE.lib;TKDEIGES.lib;TKXSBase.lib;TKXCAF.lib;TKPrim.lib;TKOpenGl.lib;TKService.lib;TKV3d.lib;TKTopAlgo.lib;TKMesh.lib;TKShHealing.lib;TKLCAF.lib;
    ```

- OCC 版本：7.9.0

  - debug版本需要自己编译，编译后文件、目录与release版本不一致

- OCC 第三方依赖：官方编译的 vc14-64 一套

  - 目前具体用到了哪些依赖可以看本目录下的`*.vcxproj.user`
    - 使用这个`*.vcxproj.user`可以实现在不同电脑中的环境独立配置，VS2022目前我只找到了这个方法
  - `CASROOT`为OCC路径，可以写死或者在环境变量中设置

- 目录：

  - `assets`: 放的是测试用的模型（后缀为.stp .step），其余的是CADFix生成的缓存
  - `include`: 放的是师兄写的几何数据结构的头文件
  - `lib`：放的是几何数据结构的库文件，（我自己重新编译过，因为与自己的NLopt库不兼容）
  - `results`：放的是debug文件（这个思路就是打印当前处理实体的类型，同步，改不改异步看实际需求）和导出的vtk、obj文件

### 项目思路

- OCC 读取模型结构有 2种（或更多）方式
  - 一种是整个模型是一个TopoCompound（混合体），然后只有一个TopoSolid
    - 不太符合导师对模型结构的要求
  - 一种是以树形结构读取，这种结构可以看CADFix、FreeCAD、Mayo（后两者都基于OCC实现，可以查看结构）。
    - 目前采用此结构，但此结构也不能满足导师的要求，因为该结构中并不是所有TopoShell都在TopoSolid中。
    - 因此我在导师推荐下，思路是 没有在TopoSolid结构下的Shell，将这些TopoShell根据**TDFLabel**分类，属于一个TDFLabel就创建一个TopoSolid
    - 所以你会发现最后**体**的数量大于CAD软件中**体**的数量
- 调用顺序
  - `importSTEP()`或`importIGES()`
  - `getBodies()`
    - `traverseLabels()`依据根标签，获取树状结构下所有TopoSolid，所有TopoShell
    - `convertToBody()`重载一 将原生TopoSolid转换为Body
    - `convertToBody()`重载二 将分类后的TopoShell，按组转换为Body
  - `convertToBody()`额，还是看代码中的注释吧，不会的再问
- 目前存在问题：
  - 记得看TODO，VS2022在最上方 视图->其他窗口->任务列表
  - `convertToNurbsSurface()`在处理圆柱体时会有问题，可能是精度问题超过了2π，然后报错
  - 目前离散采用的是一个面离散一次，不太清楚将整个模型离散，和单独将每一个面进行离散，会有哪些问题。并且也没搞清楚这个线性偏差是什么，用的时候越小，面离散的越细
  - 目前OCC自带的FixShape我没太清楚其原理，效果肯定是达不到导师要求的效果

### 资料参考

- OCC官方指南 + 官方API文档 + 官方论坛
- 大佬博客：[eryar - C++博客](http://www.cppblog.com/eryar/)（很多OCC相关中文资料出自此）
- 知乎大佬，以及他整理的资料[【opencascade】常用资料整合 - 知乎](https://zhuanlan.zhihu.com/p/448986691)
- 可以试试用Python中OCC包实现想法，再在C++中编程

