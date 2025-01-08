import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# 设置中文字体
plt.rcParams['font.sans-serif'] = ['SimHei']  # 用来正常显示中文标签
plt.rcParams['axes.unicode_minus'] = False    # 用来正常显示负号

# 模拟数据（实际使用时替换为真实数据）
# 假设片段数量在1-10之间分布
fragment_counts = np.random.normal(5, 1.5, 1000).clip(1, 10)
# 假设分子质量在100-500之间分布
molecular_masses = np.random.normal(300, 50, 1000).clip(100, 500)

# 创建图形
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

# 绘制片段数量分布
sns.histplot(data=fragment_counts, bins=20, ax=ax1, color='skyblue')
ax1.set_title('分子拆分后片段数量分布')
ax1.set_xlabel('片段数量')
ax1.set_ylabel('频数')

# 绘制分子质量分布
sns.histplot(data=molecular_masses, bins=20, ax=ax2, color='lightgreen')
ax2.set_title('分子重组后质量分布')
ax2.set_xlabel('分子质量')
ax2.set_ylabel('频数')

# 调整布局
plt.tight_layout()

# 显示图形
plt.show()

# 输出一些基本统计信息
print("\n基本统计信息：")
print("\n片段数量统计：")
print(f"平均片段数：{np.mean(fragment_counts):.2f}")
print(f"最大片段数：{np.max(fragment_counts):.2f}")
print(f"最小片段数：{np.min(fragment_counts):.2f}")

print("\n分子质量统计：")
print(f"平均分子质量：{np.mean(molecular_masses):.2f}")
print(f"最大分子质量：{np.max(molecular_masses):.2f}")
print(f"最小分子质量：{np.min(molecular_masses):.2f}")