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


library(hexSticker)
imgurl = "canOpener.svg"
sticker(imgurl,
        package = "prolfqua",
        p_color = "black",
        p_size = 7,
        p_y = 1.6,
        s_x = 1,
        s_y = 0.8,
        s_width = 0.8,
        h_fill = "green",
        h_color = "darkgreen",
        filename = "imgfile.png",
        url	= "github.com/fgcz/prolfqua",
        u_size = 2,
        u_color = "black")


#imgurl <- system.file("figures/cat.png", package="hexSticker")
sticker(imgurl, package="hexSticker", p_size=20, s_x=1, s_y=.75, s_width=.6,
        filename="cat.png")
