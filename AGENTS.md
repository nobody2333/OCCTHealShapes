# 长期工作规则

- 本文件不启动功能开发；以当前用户请求确定本轮任务。下文的具体功能需求仅在用户明确要求实施时适用。
- 咨询或研究问题不得仅据本文件修改代码；只有回答需要核实项目实现时才读取相关源码。
- 修改前检查 `git diff`，保留用户已有的未提交改动。

1. 请修改当前 C++/OCCT 项目中的 HealNarrowGap，增加对“两个 face 发生局部穿插后形成的窄边界”的修复能力。

   开始修改前，请先阅读：

   - HealNarrowGap.h / HealNarrowGap.cpp
   - GeometryUtils.h / GeometryUtils.cpp
   - TopologyUtils.h / TopologyUtils.cpp
   - RepairDiagnostics.h / RepairDiagnostics.cpp
   - ShapeRepairPipeline 相关代码
   - HealSimpleHole 和 HealClosedNarrowGap 中现有的补面、pcurve、trial shape 与校验实现
   - feature-docs/HealNarrowGap.md
- 当前 git diff
   
工作区可能已有 HealSimpleHole 和 main.cpp 的未提交修改，请保留所有无关改动，不要覆盖或格式化无关文件。
   
   # 一、问题背景
   
当前存在一种能够被 HealNarrowGap 候选检测识别的双边窄缝：
   
- 两条候选 edge 具有相同的一对拓扑端点；
   - 分别属于两个不同 face；
   - 两条 edge 距离很近；
   - 但缺陷不是普通的“两个相邻 face 边界未共享”，而是两个 face 局部发生穿插或交叉后形成的狭长边界。
   
   现有参数无法可靠修复：
   
   - ReplaceOneSide 可能因为保留 edge 不贴合另一 face 的支撑面而无法生成合法 pcurve；
- MidlineCollapse 生成的中线通常不位于任一原支撑面上；
   - SurfaceIntersection 对曲面求交精度、交线完整性和交线形态要求过高，在实际 CADFix 场景中不够稳定。
   
   本次放弃把 face-face 精确交线作为主要修复方法。
   
   # 二、总体目标
   
   复用 HealNarrowGap 当前的候选检测，不新增大型候选检测系统，也不要新建独立修复模块。
   
   在 HealNarrowGap 内增加一个“基于可信已有边界修复缺陷 face”的内部修复路径：
   
   1. 选择一条已有 edge 作为 canonical kept edge。
   2. 将另一侧 face 视为需要修复的 defective face。
   3. 用 kept edge 替换 defective face wire 中的原 edge。
   4. 优先保留 defective face 原支撑面，只修复 wire、pcurve 和必要容差。
   5. 如果原支撑面无法合法容纳修正后的边界，再根据完整修正边界重建 defective face 的支撑面。
   6. 两级都失败时拒绝候选，不强行放大容差。
   
   不要增加用户可配置的“两级修复模式”或额外开关。执行顺序固定为：
   
   原支撑面修复优先
       ↓
   仅在适合回退的失败情况下重建支撑面
   
   # 三、控制修改范围，避免过度设计
   
   本次只支持 HealNarrowGap 当前已经能够检测的场景：
   
   - 两条候选 edge 共享同一对端点；
   - 每条 edge 只属于一个 face；
   - 只替换 defective face wire 中的一条候选 edge；
   - 不实现通用 face 布尔裁剪；
   - 不计算精确 face-face 交线；
   - 不实现交线碎片拼接；
   - 不切分相邻 edge；
   - 不处理一条 edge 对多条 edge；
   - 不处理 seam、退化边、自闭合边或缺少稳定 3D curve 的边；
   - 如果 face 含有当前实现难以安全处理的复杂 inner wire，可以保守拒绝，不要扩展成通用补面框架。

   优先复用现有函数，不进行无关重构。
   
   # 四、修复入口
   
   建议增加内部函数，例如：

   bool tryRepairByRebuildingCrossedFaceBoundary(
       const TopoDS_Shape& currentShape,
       const Candidate& candidate,
       TopoDS_Shape& repairedShape,
       std::string& errorMessage);
   
   该路径内部自动执行两级修复，不向 Params 增加模式枚举或 enable 开关。
   
   可信 edge 与 defective face 的选择优先复用现有 chooseEdgeToKeep()，不要在本次建立复杂评分系统。如果选定方向无法安全修复，可以拒绝候选；不要为了提高覆盖率同时实现复杂的双向组合搜索。
   
   # 五、第一级：保留原支撑面
   
   第一级使用 canonical kept edge 替换 defective face 的原候选 edge，但继续使用 defective face 原来的 support surface。
   
   需要完成：
   
   1. 找到 kept edge 和 defective edge 在各自 wire 中的实际 edge use。
   2. 构造包含 kept edge 的修正 wire。
   3. 为 kept edge 在 defective face 原支撑面上生成 pcurve。
   4. 测量 kept edge 的 3D curve 与该支撑面/pcurve 的实际最大偏差。
   5. 只有实测偏差在独立、有限的容差预算内时才允许提高 edge tolerance。
   6. tolerance 只提高到覆盖实际误差所需的值，不得直接使用 candidateMaxGap、repairMaxGap 或两条边之间的 gap。
   7. 使用修正 wire 和原 support surface 重建 face。
   8. 保留原 face orientation、wire orientation 和非目标 wire。
9. 校验成功后提交 trial 结果。
   
   可以增加最少量、语义清晰的参数，例如：
   
   double maxBoundaryToSurfaceDeviation;
   double maxRebuiltEdgeTolerance;
   
   提供保守默认值，并在构造函数中做有限值和下限规范化。不要增加大量 filling 调参项。
   
   # 六、第一级返回状态
   
   不要通过解析 errorMessage 决定是否进入第二级。
   
   增加一个私有内部状态，例如：
   
   enum class OriginalSurfaceRepairStatus
   {
       Success,
       NeedsSurfaceRebuild,
       RejectCandidate
   };
   
   以下情况返回 NeedsSurfaceRebuild：
   
   - kept edge 到 defective face 原支撑面的实测偏差超过一级允许值；
- 无法在原支撑面上生成满足偏差要求的 pcurve；
   - 修正 wire 本身拓扑有效，但无法在原支撑面上生成有效 face。

   以下情况直接返回 RejectCandidate：
   
   - 无法可靠确定 kept edge 或 defective face；
   - 端点拓扑不一致；
   - 修正 wire 不闭合、存在分支、自交或近零长度边；
   - 存在不支持的 seam、退化边或复杂 wire；
   - 输入 shape、face、edge 或 support surface 无效；
   - 候选映射已经过期；
   - 修复将产生明显拓扑退化。

   # 七、第二级：按修正边界重建 defective face
   
   仅当第一级返回 NeedsSurfaceRebuild 时执行。
   
   第二级需要：
   
   1. 从原始 currentShape 重新创建全新的独立 trial。
   2. 重新定位 candidate、kept edge 和 defective face。
   3. 构造用 kept edge 替换 defective edge 后的完整闭合边界。
   4. 如果边界近似共面，优先创建平面 face。
   5. 非平面边界使用项目当前 OCCT 7.9.0 可用的 BRepFill_Filling、GeomPlate 或已有项目实现，以 C0 边界约束重建支撑面。
   6. 复用修正边界中的原始 edge，不为几何相同的边创建另一份独立 TShape。
   7. 为所有边界 edge 在新支撑面上建立并验证 pcurve。
   8. 检查实际边界偏差、曲面有效性和复杂度，拒绝误差过大或明显异常复杂的结果。
   9. 用新 face 替换 defective face。
   10. 保留 face 在原 shell 中的 orientation。
   
优先提取并复用 HealSimpleHole 或 HealClosedNarrowGap 中已经验证过的最小必要逻辑。不要复制整个补孔框架，不要加入渐进 degree/segments 搜索等与本缺陷无关的大型机制。采用一组保守、固定的内部构造参数即可；构造失败就拒绝候选。
   
   # 八、trial 隔离
   
   两级修复必须完全隔离：
   
- 第一级在独立 trial 上执行；
   - 第一级失败后，第二级必须从原始 currentShape 重新创建新的 trial；
   - 不得继续使用第一级失败后的 shape、face 或 edge；
   - 避免第一级对共享 TShape 的 pcurve、same-parameter 或 tolerance 修改污染第二级和原始输入。
   
   只有完整修复和最终验证通过后才能设置 repairedShape。

   # 九、SurfaceIntersection 处理

   不再让 Auto 调用 SurfaceIntersection。
   
   新的 Auto 建议顺序：

   1. 保留现有 ReplaceOneSide 快速路径；
   2. ReplaceOneSide 失败后尝试新的 crossed-face boundary repair；
   3. 如果新路径判断候选不是可安全重建的交叉 face，再按现有条件尝试 MidlineCollapse；
4. 对已经进入第二级并确认原支撑面不兼容的候选，不要再通过扩大容差或 SurfaceIntersection 强行修复。
   
   为避免扩大本次改动，可以暂时保留旧 SurfaceIntersection 枚举、参数和函数以兼容现有代码，但 Auto 不得调用它，也不要继续扩展它。若删除不会影响现有调用，可以清理；不要为了清理死代码进行大范围重构。
   
   # 十、最终拓扑要求
   
成功结果必须满足：
   
   1. kept edge 只创建或选定一次。
   2. 原正确 face 和 rebuilt defective face 共享同一个 edge TShape。
   3. 两个 face 中的 edge use 方向符合各自 wire。
   4. kept edge 在两个 face 上都具有合法 pcurve。
   5. defective edge 不再被 defective face 使用。
   6. 修正 wire 闭合且没有明显退化边。
7. rebuilt face 非空且有效。
   8. 顶层 shape 类型不变。
9. free edge 数量不增加；预期情况下应减少。
   10. 不增加新的 BRepCheck 错误类型或数量。
   11. 修复失败时最终 shape 保持不变。
   
   复用现有 validateWholeShape()、validateEdgeReplacement()、RepairDiagnostics 和 TopologyUtils；根据新 face 重建结果增加最少量必要校验。
   
   # 十一、日志
   
   增加简洁日志，至少记录：
   
   - Level 1 是否成功；
   - Level 1 失败分类：NeedsSurfaceRebuild 或 RejectCandidate；
   - 是否进入 Level 2；
   - Level 2 使用 Planar 还是 Filling；
   - 新 face 的边界偏差；
   - 最终提交或失败原因。
   
   不要逐采样点输出坐标。
   
   # 十二、文档与验证
   
   完成后：
   
   1. 更新 feature-docs/HealNarrowGap.md，说明新缺陷类型、两级自动修复流程、容差语义、限制以及删除SurfaceIntersection相关的内容。
   2. 实际编译 Debug|x64。
   3. 修复本次引入的所有编译错误和警告。
   4. 如果没有专门对应此缺陷的 STEP 模型，不新增伪造测试，也不要声称完成端到端几何验证。
   5. 可以运行现有模型做回归，但必须区分：
      - 代码修改完成；
      - 项目编译成功；
      - 现有模型回归结果；
      - 图示交叉 face 缺陷是否已经用专用模型验证。
   6. 最终报告修改文件、核心逻辑、编译结果、实际运行结果和未验证边界。
   
请直接实现完整 C++ 逻辑，不要只给设计建议或伪代码。