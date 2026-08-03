1. 请修改当前 C++/OCCT 项目中的 `HealClosedNarrowGap` 算法。

   修改目标是：

   > 保留原有对称细缝的一对一修复路径，将当前非对称细缝中“按相同归一化弧长比例建立对应关系”的方法，替换为“双向顶点投影、已有顶点吸附、单调对应节点、分段插值”的方法。
   
   当前没有提供用于验证的 STEP 文件或专门测试模型，因此本次不要求新增测试代码，也不要声称已完成端到端几何验证。修改完成后只要求项目能够实际编译成功，并保留完整的运行时安全检查。
   
   ------
   
   # 一、当前需要替换的逻辑

   当前非对称路径主要使用：

   ```cpp
   buildMidlineChainByUnifiedStations(...)
   ```

   它分别计算 chainA 和 chainB 各顶点的累计弧长比例，然后直接取两侧站点的并集，并默认：

   ```text
   chainA(t) ↔ chainB(t)
   ```
   
   即两侧相同归一化弧长位置互相对应。
   
   例如：

   ```text
   chainA: X-M-N-Z
   chainB: X-O-Z
   ```
   
   当前逻辑可能得到：
   
   ```text
   chainA stations:
   X=0.0, M=0.25, N=0.65, Z=1.0
   
   chainB stations:
   X=0.0, O=0.45, Z=1.0
   ```
   
   然后直接合并为：
   
   ```text
   0.0, 0.25, 0.45, 0.65, 1.0
   ```
   
   这种方法不能准确表达：
   
   ```text
   M 实际上对应 O
   ```
   
   但二者在各自链上的长度比例并不相同的情况。
   
   请删除或停用这种“直接合并两侧归一化顶点站点，并假设相同 t 对应”的非对称构造方式。
   
   不要删除弧长站点数据本身。弧长站点仍用于：
   
   - 表示顶点或投影点在整条链上的拓扑位置；
   - 判断投影点属于哪条 source edge；
   - 检查对应关系是否单调；
   - 构造 source edge 到 target edge 的映射。

   它不再直接决定 A、B 两侧的对应关系。
   
   ------
   
   # 二、必须保留的功能

   以下逻辑必须保留，除非为新方法进行最小必要调整：
   
   1. 自由边提取。
   2. 闭合边界环构造。
   3. 两个尖端的检测。
   4. 将边界环拆分成方向一致的 chainA 和 chainB。
   5. 原有对称细缝的一对一路径：
   
   ```cpp
   buildMidlineChainOneToOne(...)
   ```
   
   1. 原有一对一路径对应关系：
   
   ```text
   A0 ↔ B0 → M0
   A1 ↔ B1 → M1
   ```
   
   1. `EdgeReplacementTask::replacementEdgeIndices` 支持：
   
   ```text
   一条 source edge → 多条 target edges
   ```

   1. `rebuildFaceByReplacingMappedEdges()` 中：
      - 正向插入多个 replacement edges；
      - 反向时倒序并反转 replacement edges；
      - pcurve 修复；
      - wire 重建；
      - face 重建。
   2. A、B 两侧最终使用同一组 canonical target edge TShape。
   3. 内部顶点替换向 incident edges 传播的逻辑。
   4. 最终的自由边、共享边和 BRepCheck 验证。
   5. 修复失败时不修改最终 shape，只在完整修复并验证通过后提交 `trialShape`。
   
   ------
   
   # 三、修复模式命名
   
   建议将当前枚举：
   
   ```cpp
   enum class MidlineConstructionMode
   {
       LegacyOneToOne,
       UnifiedAsymmetricStations
   };
   ```

   修改为更能反映对应方式的名称，例如：

   ```cpp
   enum class MidlineConstructionMode
   {
       PairwiseEdgeMapping,
       ProjectionCorrespondence
   };
   ```
   
   含义：

   ```text
   PairwiseEdgeMapping
   原有一对一 edge 配对方式。
   
   ProjectionCorrespondence
   通过双向顶点投影建立非对称链对应关系。
   ```
   
   同步修改：
   
   ```cpp
   midlineConstructionModeName(...)
   ```
   
   以及相关诊断日志。
   
   如果为了减少改动决定保留 `LegacyOneToOne`，也可以只将：

   ```cpp
   UnifiedAsymmetricStations
   ```
   
   改名为：

   ```cpp
   ProjectionCorrespondence
   ```
   
   不要继续使用容易让人误解的 `UnifiedAsymmetricStations`。

   ------

   # 四、新增参数
   
   在 `HealClosedNarrowGap::Params` 中新增独立参数：

   ```cpp
   double projectionVertexSnapTolerance;
   ```

   该参数表示：
   
   > 投影点距离对侧已有顶点小于该值时，将投影点吸附到已有顶点。
   
   不要直接使用 `maxGapWidth` 作为顶点吸附容差，因为二者含义不同：

   ```text
   maxGapWidth
   两条细缝边链之间允许的横向距离。
   
   projectionVertexSnapTolerance
   投影点是否应当视为落在已有顶点上的容差。
   ```

   在构造函数中规范该参数：

   ```cpp
   params_.projectionVertexSnapTolerance =
       std::max(
           params_.projectionVertexSnapTolerance,
           Precision::Confusion());
   ```
   
   实际使用时限制为：
   
   ```cpp
   const double snapTolerance = std::min(
       params_.projectionVertexSnapTolerance,
       params_.maxGapWidth);
   ```
   
   如果当前项目不希望改变外部参数初始化代码，可以为该参数提供安全的默认值。
   
   不要改变现有参数的默认行为和含义。
   
   ------
   
   # 五、保留并完善链弧长数据
   
   保留或调整当前：
   
   ```cpp
   struct ChainStationData
   {
       std::vector<double> edgeLengths;
       std::vector<double> vertexStations;
       double totalLength = 0.0;
   };
   ```

   `buildNormalizedVertexStations()` 继续负责：

   1. 计算每条 source edge 的实际弧长。
   2. 计算每个 source vertex 在整条链上的累计弧长站点。
   3. 保证：
   
   ```text
   vertexStations.front() == 0.0
   vertexStations.back() == 1.0
   ```
   
   1. 保证站点严格递增。
   2. 拒绝零长度 edge、非有限长度和非法 topology count。
   
   这些站点不再直接合并成目标站点，而是用于投影结果定位和顺序判断。
   
   ------

   # 六、新增链投影结果结构

   新增结构，例如：

   ```cpp
   struct ChainProjectionResult
   {
       gp_Pnt projectedPoint;
   
       std::size_t targetEdgeIndex =
           std::numeric_limits<std::size_t>::max();
   
       Standard_Real targetParameter = 0.0;
   
       double targetStation = 0.0;
       double distance = 0.0;
   
       std::size_t snappedTargetVertexIndex =
           std::numeric_limits<std::size_t>::max();
   
       bool snappedToVertex = false;
   };
   ```

   含义：

   - `projectedPoint`：点在目标链上的最近投影点；
   - `targetEdgeIndex`：投影落在哪条 target chain edge；
   - `targetParameter`：投影点在该 edge 上的实际曲线参数；
   - `targetStation`：投影点在整条目标链上的累计弧长比例；
   - `distance`：源点到投影点的距离；
   - `snappedTargetVertexIndex`：吸附到哪个已有顶点；
   - `snappedToVertex`：是否吸附到已有顶点。
   
   所有下标必须使用明确的无效值，并在使用前检查。
   
   ------
   
   # 七、新增整条 edge chain 投影函数
   
   新增辅助函数，例如：

   ```cpp
   bool projectPointToChain(
       const gp_Pnt& sourcePoint,
       const EdgeChain& targetChain,
       const ChainStationData& targetStationData,
       ChainProjectionResult& projection,
       std::string& errorMessage) const;
   ```
   
   该函数必须将一个三维点投影到整条 `EdgeChain`，而不是只投影到某一条预先指定的 edge。
   
   实现要求：
   
   1. 遍历目标链的每条 edge。
   2. 投影范围必须限制在 edge 的实际裁剪参数范围内。
   3. 不得投影到底层无限曲线的延长部分。
   4. 正确处理 `TopLoc_Location`。
   5. 可以使用项目和当前 OCCT 版本中合适的投影接口，例如：
   
   ```text
   GeomAPI_ProjectPointOnCurve
   Extrema_ExtPC
   BRepAdaptor_Curve
   ```
   
   请根据项目当前 OCCT 版本和已有工具函数选择可编译、可靠的实现，不要强行使用不存在的接口签名。
   
   1. 对每条 edge 取得合法投影解后，选择三维距离最小的投影。
   2. 如果所有 edge 都没有有效投影，安全返回失败。
   3. 投影点距离必须满足：
   
   ```cpp
   projection.distance <=
       params_.maxGapWidth
       + params_.midlineRepairTolerance;
   ```
   
   否则当前候选不适合修复。
   
   1. 根据：
   
   ```text
   目标链中投影 edge 之前的累计长度
   +
   该 edge 起点到投影参数的实际弧长
   ```
   
   计算：
   
   ```cpp
   projection.targetStation =
       accumulatedLength / targetStationData.totalLength;
   ```
   
   1. 不要用底层曲线参数线性比例代替实际弧长比例。
   2. `targetStation` 必须限制在 `[0,1]`，且为有限数值。
   
   ------
   
   # 八、投影点吸附到已有顶点

   投影成功后，判断投影点是否靠近目标链已有顶点。

   建议只检查：

   - 投影所在 edge 的两个端点；
   - 必要时检查整个 target chain 的顶点列表。
   
   吸附条件不能只检查三维距离，还应检查沿链站点是否一致。
   
   建议逻辑：

   ```cpp
   const double stationSnapTolerance =
       std::min(
           0.25,
           snapTolerance /
               std::max(
                   targetStationData.totalLength,
                   Precision::Confusion()));
   ```

   某个目标顶点满足：

   ```cpp
   projectedPoint.Distance(targetVertexPoint)
       <= snapTolerance
   ```
   
   并且：
   
   ```cpp
   std::abs(
       projection.targetStation
       - targetStationData.vertexStations[targetVertexIndex])
       <= stationSnapTolerance
   ```

   才允许吸附。
   
   吸附后必须：
   
   ```cpp
   projection.projectedPoint =
       BRep_Tool::Pnt(targetChain.vertices[targetVertexIndex]);
   
   projection.targetStation =
       targetStationData.vertexStations[targetVertexIndex];
   
   projection.snappedTargetVertexIndex =
       targetVertexIndex;
   
   projection.snappedToVertex = true;
   ```
   
   这样能够避免把空间上很近、但在链拓扑顺序上相距很远的顶点错误吸附。
   
   ------
   
   # 九、新增对应节点结构
   
   新增非对称投影对应节点结构，例如：
   
   ```cpp
   struct ProjectionCorrespondenceNode
   {
       double stationA = 0.0;
       double stationB = 0.0;
   
       gp_Pnt pointA;
       gp_Pnt pointB;
   
       std::size_t sourceVertexIndexA =
           std::numeric_limits<std::size_t>::max();
   
       std::size_t sourceVertexIndexB =
           std::numeric_limits<std::size_t>::max();
   
       double projectionDistance = 0.0;
   
       bool hasSourceVertexA = false;
       bool hasSourceVertexB = false;
   
       bool snappedPair = false;
   };
   ```
   
   每个节点表示：
   
   ```text
   chainA(stationA) ↔ chainB(stationB)
   ```

   不再要求：

   ```text
   stationA == stationB
   ```

   目标顶点坐标最终取：
   
   ```text
   midpoint(pointA, pointB)
   ```

   ------
   
   # 十、双向投影生成对应节点
   
   新的非对称路径必须进行双向投影。
   
   ## 1. 固定端点
   
   首先创建两个强制节点：
   
   ```text
   起点：
   stationA = 0.0
   stationB = 0.0
   
   终点：
   stationA = 1.0
   stationB = 1.0
   ```

   起点和终点继续复用两侧共同的原始尖端顶点。

   ## 2. chainA 顶点投影到 chainB

   对 chainA 的每个内部顶点：

   ```cpp
   for (std::size_t vertexIndex = 1;
        vertexIndex + 1 < chainA.vertices.size();
        ++vertexIndex)
   ```

   执行：

   1. 获取源顶点坐标 `VA`。
   2. 将 `VA` 投影到整个 chainB。
   3. `stationA` 使用 A 顶点原始弧长站点。
   4. `stationB` 使用投影结果的 target station。
   5. 如果投影点吸附到 chainB 的已有顶点 `VB`：
   
   ```text
   pointA = VA
   pointB = VB
   ```
   
   并设置：
   
   ```text
   hasSourceVertexA = true
   hasSourceVertexB = true
   sourceVertexIndexA = 当前A顶点
   sourceVertexIndexB = 被吸附的B顶点
   snappedPair = true
   ```
   
   目标顶点将使用：
   
   ```text
   midpoint(VA, VB)
   ```
   
   1. 如果投影点没有吸附到已有顶点：
   
   ```text
   pointA = VA
   pointB = projectedPointOnB
   ```

   设置：
   
   ```text
   hasSourceVertexA = true
   hasSourceVertexB = false
   ```

   目标顶点将使用：

   ```text
   midpoint(VA, projectedPointOnB)
   ```
   
   ## 3. chainB 顶点投影到 chainA

   对 chainB 的每个内部顶点执行对称操作：

   ```text
   VB → project to chainA
   ```
   
   如果投影点靠近 A 的已有顶点 `VA`，建立真实顶点配对。

   否则使用：

   ```text
   midpoint(projectedPointOnA, VB)
   ```

   作为目标顶点位置。

   不能只做单向投影，因为否则可能遗漏仅存在于另一侧链上的重要内部顶点。
   
   ------
   
   # 十一、对应节点合并
   
   双向投影后可能产生重复节点。
   
   例如：
   
   ```text
   M → O
   O → M
   ```
   
   这两个投影实际上描述同一个对应关系，应合并为一个节点。
   
   新增函数，例如：

   ```cpp
   bool mergeProjectionCorrespondenceNodes(
       std::vector<ProjectionCorrespondenceNode>& nodes,
       const EdgeChain& chainA,
       const EdgeChain& chainB,
       const ChainStationData& stationDataA,
       const ChainStationData& stationDataB,
       std::string& errorMessage) const;
   ```
   
   合并条件至少包括：
   
   ```cpp
   std::abs(node1.stationA - node2.stationA)
       <= stationMergeToleranceA
   
   std::abs(node1.stationB - node2.stationB)
       <= stationMergeToleranceB
   ```

   并且目标中点距离足够近。

   站点合并容差应根据几何容差和链总长度计算，不要继续只使用：
   
   ```cpp
   128.0 * std::numeric_limits<double>::epsilon()
   ```

   建议：
   
   ```cpp
   const double stationMergeToleranceA =
       std::max(
           128.0 * std::numeric_limits<double>::epsilon(),
           snapTolerance /
               std::max(
                   stationDataA.totalLength,
                   Precision::Confusion()));
   
   const double stationMergeToleranceB =
       std::max(
           128.0 * std::numeric_limits<double>::epsilon(),
           snapTolerance /
               std::max(
                   stationDataB.totalLength,
                   Precision::Confusion()));
   ```
   
   同时设置合理上限，避免短链导致归一化容差过大。
   
   合并规则：

   1. 如果两个节点分别来自同一对真实 A、B 顶点，合并为一个 snapped pair。
   2. 如果一个是 snapped pair，另一个只是普通投影，并且二者站点接近，优先保留 snapped pair 的几何位置。
   3. 如果两个节点在同一侧包含不同的 source vertex index，不得无条件合并。
   4. 如果同一个 A 源顶点被映射到两个不同节点，返回失败。
   5. 如果同一个 B 源顶点被映射到两个不同节点，返回失败。
   6. 合并后保证每个 source internal vertex 恰好映射到一个目标节点。
   
   ------
   
   # 十二、处理多个顶点吸附到同一个对侧顶点
   
   必须处理以下冲突：
   
   ```text
   M → O
   N → O
   ```
   
   不能简单让 M、N 同时映射到同一个目标顶点，否则会导致：

   - 两个不同 A 顶点映射成同一点；
   - 目标区间长度为零；
   - 多条 source edge 被压到一个 target vertex；
   - wire 退化。
   
   处理原则：
   
   1. 一个已有目标顶点默认只能和对侧一个内部顶点形成 snapped pair。
   2. 如果多个源顶点都希望吸附到同一个目标顶点：
      - 优先保留投影距离最小的配对；
      - 另一个源顶点恢复为未吸附投影点；
      - 如果恢复后仍与已有节点过近并导致零长度区间，则拒绝候选。
   3. 不允许两个相距明显的 source vertex 最终映射到同一个 target vertex。
   4. 所有冲突处理必须确定性，不能依赖容器遍历的偶然顺序。

   ------
   
   # 十三、对应节点单调性检查
   
   所有节点合并完成后：
   
   1. 按 `stationA` 从小到大排序。
   2. 起点必须为：

   ```text
   (0.0, 0.0)
   ```
   
   1. 终点必须为：
   
   ```text
   (1.0, 1.0)
   ```
   
   1. 对任意相邻节点，必须满足：
   
   ```cpp
   next.stationA > current.stationA
   next.stationB > current.stationB
   ```

   并且差值超过对应站点容差。
   
   也就是说，对应关系必须同时沿 chainA 和 chainB 正向前进。
   
   如果出现：
   
   ```text
   M → B链70%
   N → B链40%
   ```
   
   这种投影回退，必须拒绝当前 candidate。
   
   不要尝试通过排序 stationB 或强行交换节点修复，因为这可能建立错误拓扑对应。
   
   1. 如果两个节点在任意一侧的站点过近，导致后续目标 edge 近似零长度，也应拒绝。
   2. 检查所有 A 内部 source vertex 都被映射到某个节点。
   3. 检查所有 B 内部 source vertex 都被映射到某个节点。
   
   ------
   
   # 十四、生成目标顶点
   
   根据最终单调对应节点生成 canonical target vertices。
   
   对于起点和终点：
   
   - 复用原始共同顶点；
   - 不创建新顶点；
   - 不注册顶点替换。
   
   对于每个内部节点：
   
   ```cpp
   gp_Pnt targetPoint(
       0.5 * (node.pointA.X() + node.pointB.X()),
       0.5 * (node.pointA.Y() + node.pointB.Y()),
       0.5 * (node.pointA.Z() + node.pointB.Z()));
   ```

   目标顶点容差至少覆盖：

   - `Precision::Confusion()`；
   - `params_.midlineRepairTolerance`；
   - `pointA` 到目标中点的距离；
   - `pointB` 到目标中点的距离；
   - 节点包含的 A source vertex 原始容差；
   - 节点包含的 B source vertex 原始容差；
   - source vertex 移动到目标中点的实际距离。

   为每个 A source vertex 建立：
   
   ```cpp
   targetVertexIndicesForChainA
   ```
   
   为每个 B source vertex 建立：
   
   ```cpp
   targetVertexIndicesForChainB
   ```
   
   A、B 两侧映射必须独立。
   
   如果某节点只来自 A 顶点投影到 B edge 内部，则：
   
   - A 顶点映射到该 target vertex；
   - B 侧没有 source vertex 需要替换；
   - 但该投影站点会在 B source edge 内部形成新的逻辑分段。
   
   ------
   
   # 十五、生成目标 edge
   
   对于两个相邻对应节点：
   
   ```text
   node0:
   (sA0, sB0)
   
   node1:
   (sA1, sB1)
   ```
   
   生成一条 canonical target edge。
   
   不再使用同一个全局站点 `t` 同时采样 A、B。
   
   对局部参数：
   
   ```text
   u ∈ [0,1]
   ```
   
   分别计算：
   
   ```cpp
   const double stationA =
       stationA0 + (stationA1 - stationA0) * u;
   
   const double stationB =
       stationB0 + (stationB1 - stationB0) * u;
   ```
   
   然后：
   
   ```cpp
   sampleChainPoint(chainA, stationA, pointA);
   sampleChainPoint(chainB, stationB, pointB);
   ```

   中线 guidance point 为：

   ```cpp
   midpoint(pointA, pointB)
   ```

   这样实际建立的是：

   ```text
   chainA(sA) ↔ chainB(sB)
   ```

   而不是：
   
   ```text
   chainA(t) ↔ chainB(t)
   ```
   
   采样数量根据两侧实际区间长度计算：
   
   ```cpp
   const double segmentLengthA =
       stationDataA.totalLength
       * (stationA1 - stationA0);
   
   const double segmentLengthB =
       stationDataB.totalLength
       * (stationB1 - stationB0);
   ```
   
   复用：
   
   ```cpp
   computeMidlineSampleCount(
       segmentLengthA,
       segmentLengthB);
   ```

   继续使用现有：

   ```text
   GeomAPI_Interpolate
   BRepBuilderAPI_MakeEdge
   ```
   
   每两个相邻对应节点之间只创建一条目标 edge。
   
   区间首尾 guidance point 必须使用已创建目标顶点的坐标，确保曲线端点和拓扑顶点一致。
   
   所有目标 edge 方向统一为 chain 起点到 chain 终点。

   ------

   # 十六、生成 source edge 到 target edge 的映射

   新的每条 target edge 同时对应：

   ```text
   chainA interval:
   [sA0, sA1]
   
   chainB interval:
   [sB0, sB1]
   ```
   
   使用各自区间中点：
   
   ```cpp
   const double midpointA =
       0.5 * (sA0 + sA1);
   
   const double midpointB =
       0.5 * (sB0 + sB1);
   ```
   
   分别定位：
   
   ```text
   midpointA 属于 chainA 的哪条 source edge
   midpointB 属于 chainB 的哪条 source edge
   ```

   然后填充：

   ```cpp
   replacementsForChainA[sourceEdgeIndexA]
       .push_back(targetEdgeIndex);
   
   replacementsForChainB[sourceEdgeIndexB]
       .push_back(targetEdgeIndex);
   ```
   
   例如：
   
   ```text
   chainA:
   A0
   
   chainB:
   B0, B1, B2
   ```
   
   可能得到：
   
   ```text
   A0 → {M0, M1, M2}
   
   B0 → {M0}
   B1 → {M1}
   B2 → {M2}
   ```
   
   要求：
   
   1. 每条 source edge replacement 列表非空。
   2. replacement 下标严格递增。
   3. replacement 下标连续。
   4. 每条 target edge 在 chainA 中只属于一条 source edge。
   5. 每条 target edge 在 chainB 中只属于一条 source edge。
   6. A、B 两侧都完整覆盖目标边链。
   7. 不允许同一条 wire 中多条 source edge 都映射到完全相同的 target edge。
   
   继续复用现有：
   
   ```cpp
   validateReplacementMapping(...)
   ```
   
   根据新数据结构进行必要调整。
   
   ------
   
   # 十七、替换当前非对称构造函数
   
   将当前：
   
   ```cpp
   buildMidlineChainByUnifiedStations(...)
   ```
   
   替换为类似：
   
   ```cpp
   buildMidlineChainByProjectionCorrespondence(...)
   ```
   
   建议流程：
   
   ```cpp
   bool HealClosedNarrowGap::
   buildMidlineChainByProjectionCorrespondence(
       const ClosedGapCandidate& candidate,
       MidlineChainData& midlineData,
       std::string& errorMessage) const
   {
       // 1. 构造A、B链弧长数据。
       // 2. 建立起止端点对应节点。
       // 3. A内部顶点投影到B链。
       // 4. B内部顶点投影到A链。
       // 5. 投影点吸附到已有顶点。
       // 6. 处理吸附冲突。
       // 7. 合并双向重复节点。
       // 8. 检查对应节点严格单调。
       // 9. 创建canonical target vertices。
       // 10. 创建canonical target edges。
       // 11. 构造A、B source edge映射。
       // 12. 构造A、B source vertex映射。
       // 13. 验证全部映射。
   }
   ```
   
   不要保留两套完整的非对称构造逻辑。
   
   旧的“站点并集 + 相同 t 配对”实现应删除，避免以后误用。
   
   ------
   
   # 十八、对称路径保持不变
   
   当：
   
   ```cpp
   canUseLegacyOneToOnePath(candidate)
   ```
   
   返回 true 时，继续调用现有一对一路径。
   
   不要让新增投影路径无条件接管原有对称候选。
   
   原因：
   
   - 原有对称路径已经稳定；
   - 可以保持已有输出拓扑；
   - 避免新投影方式改变原有修复结果；
   - 当前没有测试模型，不适合在本次修改中同时重构全部路径。
   
   `buildMidlineChain()` 统一入口应类似：
   
   ```cpp
   if (canUseLegacyOneToOnePath(candidate)) {
       return buildMidlineChainOneToOne(
           candidate,
           midlineData,
           errorMessage);
   }
   
   return buildMidlineChainByProjectionCorrespondence(
       candidate,
       midlineData,
       errorMessage);
   ```
   
   ------
   
   # 十九、候选检测范围
   
   本次先不要大规模修改：
   
   ```cpp
   extractCandidates(...)
   ```

   现有候选提取、尖端识别、宽度和细长度评分保持不变，以降低回归风险。

   但需要注意并在最终说明中明确：
   
   > 当前候选提取阶段仍可能使用相同归一化弧长位置进行粗略宽度判断；本次只替换非对称修复阶段的精确对应关系构造。
   
   如果 `extractCandidates()` 中存在一个与成员函数同名的局部 `sampleChainPoint` lambda，建议删除该重复实现，直接复用成员函数：
   
   ```cpp
   sampleChainPoint(...)
   ```
   
   以保证候选检测和修复阶段都使用真实弧长采样。
   
   除此之外，不要在本次修改中重写候选提取算法。
   
   ------
   
   # 二十、复用后半段重建逻辑
   
   以下函数应尽量原样复用：
   
   ```cpp
   rebuildFacesByCollapsingToMidline(...)
   rebuildFaceByReplacingMappedEdges(...)
   validateReplacementMapping(...)
   validateRepairResult(...)
   ```
   
   重点保证：
   
   1. target edges 只创建一次。
   2. A、B 两侧 face 使用相同 target edge TShape。
   3. `replacementEdgeIndices` 允许一条 source edge 插入多条 target edge。
   4. source edge 正向使用时：
   
   ```text
   M0, M1, M2
   ```
   
   1. source edge 反向使用时：
   
   ```text
   Reverse(M2), Reverse(M1), Reverse(M0)
   ```
   
   1. 每条 target edge 在每个 face 上具有正确 pcurve。
   2. `ensurePCurveOnFace()` 不得导致 A、B 两侧生成独立 target edge TShape。
   3. source internal vertex 替换继续传播到 incident edges。
   4. 所有 face 重建成功后，才将结果应用到整个 shape。

   ------
   
   # 二十一、运行时安全检查
   
   新增逻辑必须检查：
   
   1. chainA 和 chainB 非空。
   2. 两侧共享相同起点和终点。
   3. `vertices.size() == edges.size() + 1`。
   4. `ownerFaces.size() == edges.size()`。
   5. 所有 source edge 长度有效。
   6. 所有 source vertex 非空。
   7. 投影点坐标有限。
   8. 投影参数有限。
   9. 投影必须位于实际裁剪 edge 范围内。
   10. 投影距离不超过允许宽度。
   11. 投影站点在 `[0,1]` 内。
   12. source vertex 不得映射到多个不同节点。
   13. 一个对侧已有顶点不得被多个明显不同的源顶点同时吸附。
   14. 合并后的对应节点同时沿 A、B 严格单调。
   15. 不出现投影回退。
   16. 相邻节点不产生零长度 A 区间。
   17. 相邻节点不产生零长度 B 区间。
   18. target vertex 创建成功。
   19. target edge 插值成功。
   20. target edge 创建成功且非退化。
   21. source vertex 映射完整。
   22. source edge 映射完整。
   23. replacement 下标连续、递增且无重复。
   24. A、B 两侧完整覆盖 target chain。
   25. 不出现重复 target edge use。
   26. face、wire、edge、shape 均不得为空。
   27. 修复失败时不修改最终 shape。
   
   任何失败都应：
   
   ```cpp
   return false;
   ```
   
   并提供清晰、具体的 `errorMessage`。
   
   ------
   
   # 二十二、诊断日志
   
   新增简洁日志：
   
   ```text
   Midline construction mode:
   ProjectionCorrespondence
   ```
   
   投影路径至少记录：
   
   ```text
   chainAEdges
   chainBEdges
   chainAVertices
   chainBVertices
   correspondenceNodeCount
   targetVertices
   targetEdges
   ```
   
   逐个对应节点记录：
   
   ```text
   node 0:
   stationA=0.0
   stationB=0.0
   endpoint=true
   node 1:
   stationA=0.25
   stationB=0.45
   sourceVertexA=1
   sourceVertexB=1
   snappedPair=true
   projectionDistance=...
   ```
   
   或者：
   
   ```text
   node 2:
   stationA=0.65
   stationB=0.78
   sourceVertexA=2
   sourceVertexB=none
   snappedPair=false
   projectionDistance=...
   ```
   
   继续记录 edge 映射：
   
   ```text
   chainA edge 0 -> target edges [0,1,2]
   
   chainB edge 0 -> target edges [0]
   chainB edge 1 -> target edges [1]
   chainB edge 2 -> target edges [2]
   ```
   
   不要逐采样点输出坐标，避免日志过大。
   
   ------
   
   # 二十三、中文注释要求
   
   在新增核心代码处增加简洁、准确的中文注释，重点说明：
   
   1. 为什么投影范围必须限制在裁剪 edge 内。
   2. 为什么投影点还要转换成累计弧长站点。
   3. 为什么顶点吸附容差不能直接等同于最大细缝宽度。
   4. 为什么需要双向投影。
   5. 为什么 M 投影到 O 附近时，应使用 M 和 O 的中点。
   6. 为什么未吸附时使用源顶点与投影点的中点。
   7. 为什么多个源顶点不能同时吸附到同一个对侧顶点。
   8. 为什么对应节点必须同时沿 A、B 单调递增。
   9. 为什么目标 edge 内部采样需要分别插值 `stationA` 和 `stationB`。
   10. 为什么 target edge 只创建一次并由两侧 face 共享。
   11. 为什么一条 source edge 可以映射到多条 target edge。
   12. 为什么反向 source edge 需要倒序并反转 replacement edges。
   
   不要为普通赋值和显而易见的循环添加无意义注释。
   
   ------
   
   # 二十四、文件修改范围
   
   请先阅读：
   
   ```text
   HealClosedNarrowGap.h
   HealClosedNarrowGap.cpp
   GeometryUtils.h
   GeometryUtils.cpp
   TopologyUtils.h
   TopologyUtils.cpp
   RepairDiagnostics.h
   RepairDiagnostics.cpp
   相关调用代码
   Visual Studio项目配置
   ```
   
   预计主要修改：
   
   ```text
   HealClosedNarrowGap.h
   HealClosedNarrowGap.cpp
   ```
   
   如果投影或弧长计算适合放入 `GeometryUtils`，可以进行最小必要修改。
   
   不要对无关文件做大范围格式化或重构。
   
   ------
   
   # 二十五、编译要求
   
   当前没有测试文件，因此本次不要求新增或运行几何测试。
   
   修改完成后必须实际编译项目。
   
   要求：
   
   1. 修复本次修改引入的所有编译错误。
   2. 修复本次修改引入的所有编译警告。
   3. 项目可能启用了 `/WX`，因此注意：
      - `int` 与 `std::size_t` 比较；
      - 有符号和无符号转换；
      - 未使用变量；
      - 未初始化下标；
      - lambda 返回类型；
      - 缺少必要头文件；
      - OCCT接口签名；
      - `Standard_Real` 与 `double` 转换；
      - `TopLoc_Location` 处理；
      - `Handle(Geom_Curve)` 的空值检查；
      - 无效枚举分支；
      - 重复函数定义。
   4. 不得在没有实际编译时声称编译通过。
   5. 如果由于环境、依赖库或项目配置无法编译，明确说明阻塞原因。
   6. 不要新增依赖外部 STEP 文件的测试代码。
   7. 不要声称新投影算法已经完成端到端几何验证。
   8. 明确区分：
   
   ```text
   代码修改完成
   项目编译成功
   实际几何模型验证完成
   ```
   
   本次没有测试模型，因此第三项应说明尚未完成。
   
   ------
   
   # 二十六、完成后输出
   
   修改完成后，请结合以下内容重新整理feature-docs中HealClosedNarrowGap.md文件：
   
   1. 修改了哪些文件。
   2. 每个文件的核心修改内容。
   3. 删除或停用了哪些旧的非对称站点并集逻辑。
   4. 新增了哪些数据结构。
   5. 新增了哪些投影辅助函数。
   6. 点到整条 edge chain 的投影如何实现。
   7. 如何保证投影限制在实际裁剪 edge 内。
   8. 投影点如何转换成累计弧长站点。
   9. 顶点吸附容差如何使用。
   10. 双向投影如何合并。
   11. 多个源顶点争用同一对侧顶点时如何处理。
   12. 如何检查对应节点单调性。
   13. target vertex 如何生成。
   14. target edge 如何使用 `(stationA, stationB)` 分段插值生成。
   15. source edge 到 target edge 的一对多映射如何生成。
   16. 原有对称一对一路径如何保持不变。
   17. 后半段 wire、pcurve 和 face 重建如何复用。
   18. 增加了哪些中文注释和运行时安全检查。
   19. 实际编译结果。
   20. 当前尚未验证或尚未覆盖的边界情况。
   
   不要只给修改建议，请直接修改代码。
   
   不要使用伪代码代替以下核心实现：
   
   - 点到整条链的投影；
   - 投影点弧长站点计算；
   - 顶点吸附；
   - 双向对应节点生成；
   - 节点合并；
   - 冲突处理；
   - 单调性检查；
   - target vertex 生成；
   - target edge 生成；
   - source edge 映射；
   - source vertex 映射。
   
   由于没有测试模型，不要虚构测试结果。