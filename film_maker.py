#ffmpeg -f image2 -s 1920x1080 -start_number 10 -i dens_x_3_%05d.png -r 10 -pix_fmt yuv420p -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2:color=white" -c:v libx264 -crf 15 movie.mp4



#ffmpeg -f image2 -s 1920x1080 -start_number 10 -r 8 -i dens_x_3_%05d.png -pix_fmt yuv420p -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2:color=white" -c:v libx264 -b:v 4M movie.mp4




#Densite suivant x
ffmpeg -f image2 -s 1920x1080 -start_number 10 -r 8 -i dens_x_3_%05d.png -pix_fmt yuv420p -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2:color=white" -c:v libx264 -crf 0  movie_dens_x.mp4

#Densite suivant y
ffmpeg -f image2 -s 1920x1080 -start_number 10 -r 8 -i dens_y_3_%05d.png -pix_fmt yuv420p -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2:color=white" -c:v libx264 -crf 0  movie_dens_y.mp4

#Densite suivant z
ffmpeg -f image2 -s 1920x1080 -start_number 10 -r 8 -i dens_z_3_%05d.png -pix_fmt yuv420p -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2:color=white" -c:v libx264 -crf 0  movie_dens_z.mp4


#Vitesse suivant x
ffmpeg -f image2 -s 1920x1080 -start_number 10 -r 8 -i vel_x_3_%05d.png -pix_fmt yuv420p -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2:color=white" -c:v libx264 -crf 0  movie_vel_x.mp4

#Vitesse suivant y
ffmpeg -f image2 -s 1920x1080 -start_number 10 -r 8 -i vel_y_3_%05d.png -pix_fmt yuv420p -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2:color=white" -c:v libx264 -crf 0  movie_vel_y.mp4

#Vitesse suivant z
ffmpeg -f image2 -s 1920x1080 -start_number 10 -r 8 -i vel_z_3_%05d.png -pix_fmt yuv420p -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2:color=white" -c:v libx264 -crf 0  movie_vel_z.mp4
