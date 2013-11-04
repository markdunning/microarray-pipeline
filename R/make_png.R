make_png <-
function (file_name){
###Just a quick function to make pngs, they you can globally set the image size and options
    png(file_name,width=as.numeric(config["png_w",]),height=as.numeric(config["png_h",]),type = "cairo")
  }
