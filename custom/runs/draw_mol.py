import io
from rdkit import Chem
from rdkit.Chem import Draw
from PIL import Image, ImageDraw, ImageFont


from rdkit import Chem
from rdkit.Chem import Draw
from PIL import Image, ImageDraw, ImageFont
import io


def draw_molecule(smiles, img_size=(600, 600), dpi=400):
    # 将SMILES字符串转换为分子对象
    molecule = Chem.MolFromSmiles(smiles)

    # 绘制分子图并保存为PIL图像
    drawer = Draw.MolDraw2DCairo(img_size[0], img_size[1])
    drawer.drawOptions().fixedBondLength = 40  # 调整键的长度
    drawer.DrawMolecule(molecule)
    drawer.FinishDrawing()
    mol_img = Image.open(io.BytesIO(drawer.GetDrawingText()))

    # 创建一个新的图像，尺寸比分子图像大，以便在底部添加文本
    width, height = mol_img.size
    total_height = height + 50  # 留出空间放置文本
    new_img = Image.new("RGB", (width, total_height), "white")

    # 将分子图粘贴到新图像的顶部
    new_img.paste(mol_img, (0, 0))

    # 在底部添加SMILES序列文本
    draw = ImageDraw.Draw(new_img)
    try:
        font = ImageFont.truetype("arial.ttf", 30)  # 使用指定字体文件和大小
    except IOError:
        font = ImageFont.load_default()  # 回退到默认字体
    
    try:
        text_width, text_height = draw.textbbox((0, 0), smiles, font=font)[2:]
    except AttributeError:
        text_width, text_height = draw.textsize(smiles, font=font)

    text_position = ((width - text_width) // 2, height + 5)
    draw.text(text_position, smiles, fill="black", font=font)

    return new_img

def create_molecule_grid(smiles_list, grid_size=(2, 2), border_width=5, border_color="black"):
    images = [draw_molecule(smiles) for smiles in smiles_list]

    # 获取单个图像的尺寸
    single_width, single_height = images[0].size
    grid_width = single_width * grid_size[0]
    grid_height = single_height * grid_size[1]

    # 创建一个新的图像来包含所有的分子图像
    grid_img = Image.new("RGB", (grid_width, grid_height), "white")
    draw = ImageDraw.Draw(grid_img)

    # 将每个分子图像粘贴到网格中的适当位置并绘制边框
    for i, img in enumerate(images):
        x = (i % grid_size[0]) * single_width
        y = (i // grid_size[0]) * single_height
        grid_img.paste(img, (x, y))
        draw.rectangle(
            [x, y, x + single_width - 1, y + single_height - 1],
            outline=border_color,
            width=border_width
        )

    # 显示最终图像
    grid_img.show()

    # 保存图像为文件
    grid_img.save("molecule_grid.png")


# 示例SMILES字符串列表
smiles_list = [
    "[14*]c1nc([14*])c([14*])[nH]1",
    "[8*]C(F)(F)F",
    "[16*]c1ccc2ncsc2c1",
    "[14*]c1nc(C)cs1",
]

# 创建并显示分子组图
create_molecule_grid(smiles_list)
