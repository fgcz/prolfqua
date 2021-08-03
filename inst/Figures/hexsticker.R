library(hexSticker)
imgurl = "C:\\Users\\wolski\\__checkout\\prolfqua\\inst\\Figures\\Nalweka_2.jpg"
sticker(imgurl,
        package = "prolfqua",
        p_color = "yellow",
        p_size = 25,
        p_y = 1,
        s_x = 1,
        s_y = 1,
        s_width = 0.3,
        h_fill = "green", h_color = "darkgreen",
        filename = "imgfile.png",
        url	= "github.com/wolski/prolfqua",
        u_size = 5.2,
        u_color = "white")


imgurl <- system.file("figures/cat.png", package="hexSticker")
sticker(imgurl, package="hexSticker", p_size=20, s_x=1, s_y=.75, s_width=.6,
        filename="cat.png")
