# 眼科科研助手 - 永久记忆库

## 📋 个人信息与研究方向

| 项目 | 内容 |
|------|------|
| **研究者** | 小cui |
| **研究方向** | 视网膜疾病（全方向） |
| **涵盖领域** | AMD、糖尿病视网膜病变、RVO、视网膜脱离、遗传性视网膜病变 |
| **知识库初始化日期** | 2026-03-05 |

## 🎯 当前研究项目

### [项目1] - （待命名）
**状态:** 策划阶段 | **优先级:** 🔴 高
- 研究问题: （待补充）
- 疾病方向: AMD / DME / RVO / 其他
- 研究类型: （待补充）
- 目标期刊: （待补充）

### [项目2] - （待命名）
**状态:** 待启动 | **优先级:** 🟡 中

### [项目3] - （待命名）
**状态:** 待启动 | **优先级:** 🟢 低

## 📊 统计学方法论库

### Logistic回归
适用场景: 二分类结局
R实现: glm(outcome ~ var1, family = binomial)

### 线性混合效应模型
适用场景: 重复测量、双眼数据
R实现: lmer(outcome ~ var1 + (1|patient_id))

### 双眼数据处理方案
1. 单眼分析（简单）
2. 混合效应模型（推荐）
3. GEE（大样本）

## 📚 Target Journals

**Tier 1**: Ophthalmology, JAMA Ophthalmology
**Tier 2**: Retina, IOVS, British Journal of Ophthalmology
**Tier 3**: Ophthalmology Retina, Eye

## 📅 更新日志

| 日期 | 内容 | 状态 |
|------|------|------|
| 2026-03-05 | 初始化 | ✅ |
| 2026-03-07 | 文献检索功能优化，成功获取Ophthalmology和IOVS最新期刊 | ✅ |
| 2026-03-07 | 新增文献检索快捷方式 (r.jina.ai代理) | ✅ |

*最后更新: 2026-03-07*

---

## 🔧 工具配置

### 文献检索
- **主工具**: web_fetch (通过r.jina.ai代理)
- **Ophthalmology**: `https://r.jina.ai/http://www.aaojournal.org/current`
- **IOVS**: `https://r.jina.ai/http://iovs.arvojournals.org`
- **PubMed**: 当前不可用 (地区限制 EO 14117)

### 替代检索策略
1. 使用Bing/Google搜索 + jina.ai摘要服务
2. 直接访问期刊网站目录页
3. 使用Google Scholar (如可用)
