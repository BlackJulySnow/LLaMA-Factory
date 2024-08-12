import json
import matplotlib.pyplot as plt

train_losses = []
eval_losses = []
steps = []

with open("log_file.txt", "r") as file:
    for line in file:
        # 解析每一行的JSON
        data = json.loads(line.strip())
        train_losses.append(data['loss'])  # 假设文件中提供的 loss 是 train loss
        steps.append(data['current_steps'])
        
        # 如果文件中有 eval loss，可以按类似方式提取
        # 假设 eval loss 在文件中对应的 key 是 'eval_loss'
        if 'eval_loss' in data:
            eval_losses.append(data['eval_loss'])
        else:
            eval_losses.append(None)

# 绘制 train loss
plt.plot(steps, train_losses, label='Train Loss')

# 绘制 eval loss，如果存在的话
if any(eval_losses):  # 检查是否存在非 None 的 eval loss
    plt.plot(steps, eval_losses, label='Eval Loss')

plt.xlabel('Steps')
plt.ylabel('Loss')
plt.title('Train and Eval Loss over Steps')
plt.legend()
plt.show()
