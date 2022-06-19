# Generating the QRCODE
library("qrcode")
png("qrplot_project.png")
qrcode_gen("https://evasyutyicheng.shinyapps.io/Project/")
dev.off()

# Show the QRCode
library("imager")
im<-load.image("qrplot_project.png")
plot(im)