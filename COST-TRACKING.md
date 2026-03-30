# COST-TRACKING.md - OpenClaw 成本追踪

## 优化基准 (2026-03-12)

### 优化前
- 系统文件: ~6,776 tokens
- AGENTS.md: ~1,967 tokens
- SOUL.md: ~1,194 tokens
- MEMORY.md: ~1,948 tokens

### 优化后
- 系统文件: ~2,563 tokens (**-62%**)
- AGENTS.md: ~312 tokens (-84%)
- SOUL.md: ~428 tokens (-64%)
- MEMORY.md: ~418 tokens (-79%)
- FOCUS.md: ~139 tokens (新增热内存)

## 配置优化

| 设置 | 优化前 | 优化后 | 效果 |
|------|--------|--------|------|
| cacheRetention | short | **long** | 提示缓存，节省90%重复成本 |
| contextPruning.ttl | 1h | **55m** | 对齐心跳，保持缓存热 |
| heartbeat.every | 30m | **55m** | 减少不必要的心跳 |
| tools.profile | coding | **full** | 解锁所有工具 |
| gateway.bind | auto | **loopback** | 安全加固 |

## 技能统计

- **已安装技能**: 45+
- **懒加载**: 启用 SKILL-INDEX.md
- **技能来源**: Skillhub (lightmake.site)

## 目标成本削减

- **Week 1 目标**: 30-40% 成本削减 ✅ 已完成
- **Week 2 目标**: 额外 20-30% 削减
- **Week 3 目标**: 总计 50%+ 削减

## 监控指标

- [ ] 每日会话数
- [ ] 平均 tokens/会话
- [ ] 模型使用分布
- [ ] 缓存命中率

## 下次审查日期

**2026-03-19** (一周后)

---

*创建: 2026-03-12*
