options(stringsAsFactors=F)
if (interactive()) {
  setwd('~/d/sci/src/aso_in_vitro')
}
library(plotrix)
library(sqldf)
library(minpack.lm)
library(drc)
library(magick)

imgmode = 'png'
# imgmode = 'pdf'
imgsave = get(imgmode) # get the function that saves as that type of image
if (imgmode == 'png') {
  resx = 600 # multiplier for image width, height, and resolution
} else if (imgmode == 'pdf') {
  resx = 1
}

percent = function(proportion,digits=2) {
  return ( gsub(' ','',paste(formatC(proportion*100, digits=digits, format='fg'),"%",sep="") ) )
}

alpha = function(rgb_hexcolor, proportion) {
  hex_proportion = sprintf("%02x",round(proportion*255))
  rgba = paste(rgb_hexcolor,hex_proportion,sep='')
  return (rgba)
}

#### FIGURE 1: ITC
imgsave(paste('figures/figure-1.',imgmode,sep=''),width=6.5*resx,height=6*resx,res=resx)

layout_matrix = matrix(c(1,1,1,1,1,3,3,3,3,3,2,2,2,2,2,4,4,4,4,4,5,5,6,6,7,7,8,8,9,9,5,5,6,6,7,7,8,8,9,9), nrow=4, byrow=T)
layout(layout_matrix)

top_ylim = c(-1, .05)
bot_ylim = c(-1.25e5, 2e4)
xats = (0:25)/10
xbigs = (0:10)*.25

time_xlim = c(0,4863)/2
molar_xlim = c(0, 2.1125)/2

# active ASO 1
par(mar=c(0,5,3,2))
isotherm = read.table('data/itc/active_aso_1_prp1_raw.txt',sep='\t',header=T)
plot(isotherm$time_x, isotherm$cp_y, xlim=time_xlim, ylim=top_ylim, type='l', axes=F, ann=F, xaxs='i', yaxs='i')
axis(side=2, at=c(-1, -.75, -.5, -.25, 0, .25), lwd=1, lwd.ticks=1, las=2)
mtext(side=2, line=3, text='µcal/s')
mtext('A', side=3, cex=2, adj = -0.1, line = 0.5)
par(mar=c(4,5,1,2))
fits = read.table('data/itc/active_aso_1_fits.tsv',sep='\t',header=T)
colnames(fits) = c('dh','injv','xt','mt','xmt','ndh','fit','dy')
fits$ndh = as.numeric(fits$ndh)
fits$xmt = as.numeric(fits$xmt)
xlim = range(fits$xmt, na.rm=T)
plot(NA, NA, xlim=molar_xlim, ylim=bot_ylim, axes=F, ann=F, xaxs='i', yaxs='i')
abline(h=0, lwd=.25)
points(fits$xmt, fits$ndh, pch=15, cex=1, col='black')
points(fits$xmt, fits$fit,type='l',lwd=1,col='red')
axis(side=1, at=xats, labels=NA, tck=-0.025, lwd=1, lwd.ticks=1)
axis(side=1, at=xbigs, tck=-0.05, lwd=0, lwd.ticks=1)
mtext(side=1, line=2.5, text='molar ratio')
axis(side=2, at=c(-1.5e5, -1e5, -5e4, 0, 5e4), labels=c('-150','-100', '-50', '0', ''), las=2)
mtext(side=2, line=3, text='kcal/mol')

mass_time_xlim = time_xlim*1.2
mass_xlim = molar_xlim*1.2
par(mar=c(0,5,3,2))
isotherm = read.table('data/itc/heparin_prp_9raw.txt',sep='\t',header=T)
plot(isotherm$time_x, isotherm$cp_y, xlim=mass_time_xlim, ylim=top_ylim, type='l', axes=F, ann=F, xaxs='i', yaxs='i')
axis(side=2, at=c(-1, -.75, -.5, -.25, 0, .25), lwd=1, lwd.ticks=1, las=2)
mtext(side=2, line=3, text='µcal/s')
mtext('B', side=3, cex=2, adj = -0.1, line = 0.5)
par(mar=c(4,5,1,2))
fits = read.table('data/itc/heparin_fits.tsv',sep='\t',header=T)
colnames(fits) = c('dh','injv','xt','mt','xmt','ndh','fit','dy')
# FYI: in these files, ndh is background-subtracted kcal/mol, and xmt is molar ratio
fits$ndh = as.numeric(fits$ndh)
fits$xmt = as.numeric(fits$xmt)
xlim = range(fits$xmt, na.rm=T)
plot(NA, NA, xlim=mass_xlim, ylim=bot_ylim, axes=F, ann=F, xaxs='i', yaxs='i')
abline(h=0, lwd=.25)
points(fits$xmt, fits$ndh, pch=15, cex=1, col='black')
points(fits$xmt, fits$fit,type='l',lwd=1,col='red')
axis(side=2, at=c(-1.5e5, -1e5, -5e4, 0, 5e4), labels=c('-150','-100', '-50', '0', ''), las=2)
mtext(side=2, line=3, text='kcal/mol')
# for heparin it is hetereogeneous molecular weight so convert from the pseudo-molar ratio
# based on 10kDa into a mass ratio
# to convert back from molar ratio use PrP: 22.7 kDa, hep: 10 kDa
# mass ratio = molar ratio * 10/22.7
mass_conversion = 10/22.7
fits$massratio = fits$xmt * mass_conversion
mass_xats = (0:10)/10
mass_xbigs = (0:4)*.25
axis(side=1, at=mass_xats/mass_conversion, labels=NA, tck=-0.025, lwd=1, lwd.ticks=1)
axis(side=1, at=mass_xbigs/mass_conversion, labels=mass_xbigs, tck=-0.05, lwd=0, lwd.ticks=1)
mtext(side=1, line=2.5, text='mass ratio')

par(mar=c(4,5,3,2))
# 2 empty plots to reserve space for left hand material
plot(NA, NA, xlim=c(0,1), ylim=c(0,1), ann=F, axes=F)
plot(NA, NA, xlim=c(0,1), ylim=c(0,1), ann=F, axes=F)

# thermodynamic parameters
params = read.table('data/itc/itc_params.tsv',sep='\t',header=T)
itc = read.table('data/itc/itc_all_asos.tsv',sep='\t',header=T)
itc$kd = 1/itc$ka
itcsum = sqldf("
               select   p.y, p.display, avg(i.kd) kd_mean, stdev(i.kd) kd_sd, count(*) n,
               avg(i.h) h_mean, stdev(i.h) h_sd, avg(i.s) s_mean, stdev(i.s) s_sd
               from     itc i, params p
               where    i.aso = p.display
               group by 1
               order by 1 desc
               ;")
itcsum$kd_se = itcsum$kd_sd / sqrt(itcsum$n)
itcsum$kd_l95 = itcsum$kd_mean - 1.96 * itcsum$kd_se
itcsum$kd_u95 = itcsum$kd_mean + 1.96 * itcsum$kd_se

itcsum$h_se = itcsum$h_sd / sqrt(itcsum$n)
itcsum$h_l95 = itcsum$h_mean - 1.96 * itcsum$h_se
itcsum$h_u95 = itcsum$h_mean + 1.96 * itcsum$h_se

itcsum$s_se = itcsum$s_sd / sqrt(itcsum$n)
itcsum$s_l95 = itcsum$s_mean - 1.96 * itcsum$s_se
itcsum$s_u95 = itcsum$s_mean + 1.96 * itcsum$s_se

asocats = sqldf("
select   type, avg(y) meany, max(y) + 0.35 topy, min(y) - 0.35 boty
from     params
group by 1
order by 1
;")
asocats$label = gsub(', ','\n',asocats$type)

xlim = c(-9, -6)
ylim = c(0.5, 8.5)
xbigs = -9:-6
xbigs_labs1 = c('1 nM','10 nM','100 nM','1 µM')
xats = log10(rep(1:9,4) * rep(10^(-9:-6),each=9))

par(mar=c(4,1,3,2))
plot(NA, NA, xlim=xlim, ylim=ylim, xaxs='i', yaxs='i', ann=F, axes=F)
axis(side=1, at=xlim, labels=NA, lwd=1, lwd.ticks=0)
axis(side=1, at=xats, labels=NA, lwd=0, lwd.ticks=1, tck=-0.025)
axis(side=1, at=xbigs, labels=xbigs_labs1, lwd=0, lwd.ticks=1, tck=-0.05)
mtext(side=1, line=2.5, text='Kd')
segments(x0=log10(itcsum$kd_l95), x1=log10(itcsum$kd_u95), y0=itcsum$y)
points(x=log10(itcsum$kd_mean), y=itcsum$y,pch=19)
mtext('C', side=3, cex=2, adj = -0.1, line = 0.5)

# use the extra space to the left
par(xpd=T)
axis(side=2, at=ylim, labels=NA, lwd=1, lwd.ticks=0)
axis(side=2, at=itcsum$y, labels=itcsum$display, las=2, lwd=0, lwd.ticks=0)
for (i in 1:nrow(asocats)) {
  axis(side=2, line=8, at=c(asocats$boty[i], asocats$topy[i]), labels=NA, lwd=1, lwd.ticks=1, las=2, tck=0.05)
  axis(side=2, line=8, at=asocats$meany[i], labels=asocats$label[i], lwd=0, lwd.ticks=0, las=2, cex.axis=1.2)
}
par(xpd=F)

par(mar=c(4,1,3,2))
xlim=c(-100000,0)
xats = 10000 * (-10:0)
xbigs = c(-100000, -50000, 0)
xbigs_labs1 = c(-100, -50, 0)
plot(NA, NA, xlim=xlim, ylim=ylim, xaxs='i', yaxs='i', ann=F, axes=F)
axis(side=2, at=ylim, labels=NA, lwd=1, lwd.ticks=0)
axis(side=1, at=xlim, labels=NA, lwd=1, lwd.ticks=0)
axis(side=1, at=xats, labels=NA, lwd=0, lwd.ticks=1, tck=-0.025)
axis(side=1, at=xbigs, labels=xbigs_labs1, lwd=0, lwd.ticks=1, tck=-0.05)
mtext(side=1, line=2.5, text='ΔH (kcal/mol)')
segments(x0=itcsum$h_l95, x1=itcsum$h_u95, y0=itcsum$y)
points(x=itcsum$h_mean, y=itcsum$y,pch=19)
mtext('D', side=3, cex=2, adj = -0.1, line = 0.5)

par(mar=c(4,1,3,2))
xlim=c(-300,0)
xats = 10 * (-30:0)
xbigs = 100* (-3:0)
xbigs_labs1 = 100* (-3:0)
plot(NA, NA, xlim=xlim, ylim=ylim, xaxs='i', yaxs='i', ann=F, axes=F)
axis(side=2, at=ylim, labels=NA, lwd=1, lwd.ticks=0)
axis(side=1, at=xlim, labels=NA, lwd=1, lwd.ticks=0)
axis(side=1, at=xats, labels=NA, lwd=0, lwd.ticks=1, tck=-0.025)
axis(side=1, at=xbigs, labels=xbigs_labs1, lwd=0, lwd.ticks=1, tck=-0.05)
mtext(side=1, line=2.5, text='ΔS (cal/mol/°K)')
segments(x0=itcsum$s_l95, x1=itcsum$s_u95, y0=itcsum$y)
points(x=itcsum$s_mean, y=itcsum$y,pch=19)
mtext('E', side=3, cex=2, adj = -0.1, line = 0.5)

dev.off() #### -- END FIGURE 1












#### FIGURE 2: AGGREGATION BY NMR, DLS
imgsave(paste('figures/figure-2.',imgmode,sep=''),width=6.5*resx,height=5*resx,res=resx)

layout_matrix = matrix(c(1,1,2,3,3,3),nrow=2,byrow=T)
layout(layout_matrix)

asocol = '#EEAA00'
prpcol = '#0000CC'
mixcol = '#AA7799'

# Panel A: NMR

path = 'data/nmr/PrP8_50uM_control.matrix'
spec = as.matrix(read.table(path,row.names=1,header=TRUE,quote='',comment.char='',skip=8,sep='\t'))
spec = pmax(spec, 0.0)
spec = pmin(spec, 20000)
spec = t(spec)
control_spec = spec

path = 'data/nmr/PrP8_50uM_activeASO1_50uM.matrix'
spec = as.matrix(read.table(path,row.names=1,header=TRUE,quote='',comment.char='',skip=8,sep='\t'))
spec = pmax(spec, 0.0)
spec = pmin(spec, 20000)
spec = t(spec)
aso_spec = spec

n15vals = -as.numeric(gsub('y','',colnames(spec)))
h1vals = -as.numeric(gsub('x','',rownames(spec)))

line_color = '#777777'
axis_color = '#000000'

par(mar=c(4,6,3,2))
plot(NA, NA, axes=F, ann=F, xlim=c(-10.1, -5.9), ylim=c(-131, -106), xaxs='i', yaxs='i')
xats = -11:-4
yats = (-10:-13)*10
abline(h=yats, col=line_color,lwd=.25)
abline(v=xats, col=line_color,lwd=.25)
axis(side=1, at=xats, lwd=0, lwd.ticks=0 ,col.axis=axis_color)
axis(side=2, at=yats, lwd=0, lwd.ticks=0, las=2, col.axis=axis_color)
levels = (4:9)*1000
contour(x=h1vals, y=n15vals, control_spec, levels=levels, axes=FALSE, ann = FALSE, add=TRUE, col=prpcol)
contour(x=h1vals, y=n15vals, aso_spec, levels=levels, axes=FALSE, ann = FALSE, add=TRUE, col=mixcol)
mtext(side=1, line=2.5, text=expression(''^'1'~'H (ppm)'))
mtext(side=2, line=3.5, text=expression(''^'15'~'N (ppm)'))  
legend(x=-7,y=-105,legend=c('PrP alone','ASO + PrP'),col=c(prpcol,mixcol),pch=15,text.col=c(prpcol,mixcol),bty='n')
mtext('A', side=3, cex=2, adj = -0.1, line = 0.5)

# Panel B: using magick to display photo of NMR tubes
img = image_read('data/nmr/IMG_20170608_103051.jpg')
img_png = image_convert(img, 'png')
plot(as.raster(img_png))
mtext('B', side=3, cex=2, adj = -0.1, line = 0.5)

# Panel C: DLS histos

dls = read.table('data/dls/190425_DLS_histograms_PrP-ASO.txt',sep='\t',header=F,comment.char='',quote='')
colnames(dls) = c('radius','aso50','prp50','asoprp50','blank','rad2','prp10','asoprp10_1','asoprp10_o')
dls$aso50 = dls$aso50/100
dls$prp50 = dls$prp50/100
dls$asoprp50 = dls$asoprp50/100

ylim = c(0,.75)
xats = c((1:9)/10, 1:9, (1:9)*10, (1:9)*100, (1:9)*1000)
xbigs = c(.1, 1, 10, 100)
xbigs_labs = c('0.1','1','10','100')

plot(NA, NA, xlim=c(0.3,300), ylim=ylim, log='x', ann=F, axes=F, xaxs='i', yaxs='i')
points(dls$radius, dls$prp50, col=prpcol, type='h',lwd=20, lend=1)
points(dls$radius, dls$aso50, col=asocol, type='h',lwd=20, lend=1)
points(dls$radius, dls$asoprp50, col=mixcol, type='h',lwd=20, lend=1)
axis(side=1, at=xats, labels=NA, tck=-0.025)
axis(side=1, at=xbigs, labels=xbigs_labs, lwd=0, lwd.ticks=1, tck=-0.05)
axis(side=2, at=(0:10)/10, labels=percent((0:10)/10), las=2)
mtext(side=1, line=2.5, text='hydrodynamic radius (nm)')
mtext(side=2, line=3, text='proportion of mass')
legend('topright',legend=c('ASO alone','PrP alone','ASO + PrP'),col=c(asocol,prpcol,mixcol),pch=15,text.col=c(asocol,prpcol,mixcol),bty='n')
mtext('C', side=3, cex=2, adj = -0.1, line = 0.5)

dev.off() #### --- END FIGURE 2






#### FIGURE 3: ITC SALT DEPENDENCE
imgsave(paste('figures/figure-3.',imgmode,sep=''),width=8.5*resx,height=4*resx,res=resx)

par(mfrow=c(1,2), mar=c(4,4,3,2))

ninj = 40 # injections

a1 = read.table('data/itc-nacl/A1.txt',sep='\t')
colnames(a1) = c('dh','injv','xt','mt','xmt','ndh','fit','dy')
a1 = a1[2:ninj,]
a1$ndh = as.numeric(a1$ndh)
a1$xmt = as.numeric(a1$xmt)

a3 = read.table('data/itc-nacl/A3.txt',sep='\t')
colnames(a3) = c('dh','injv','xt','mt','xmt','ndh','fit','dy')
a3 = a3[2:ninj,]
a3$ndh = as.numeric(a3$ndh)
a3$xmt = as.numeric(a3$xmt)

a5 = read.table('data/itc-nacl/A5.txt',sep='\t')
colnames(a5) = c('dh','injv','xt','mt','xmt','ndh','fit','dy')
a5 = a5[2:ninj,]
a5$ndh = as.numeric(a5$ndh)
a5$xmt = as.numeric(a5$xmt)

a7 = read.table('data/itc-nacl/A7.txt',sep='\t')
colnames(a7) = c('dh','injv','xt','mt','xmt','ndh')
a7 = a7[2:ninj,]
a7$ndh = as.numeric(a7$ndh)
a7$xmt = as.numeric(a7$xmt)

b1 = read.table('data/itc-nacl/B1.txt',sep='\t')
colnames(b1) = c('dh','injv','xt','mt','xmt','ndh','fit','dy')
b1 = b1[2:ninj,]
b1$ndh = as.numeric(b1$ndh)
b1$xmt = as.numeric(b1$xmt)

b3 = read.table('data/itc-nacl/B3.txt',sep='\t')
colnames(b3) = c('dh','injv','xt','mt','xmt','ndh','fit','dy')
b3 = b3[2:ninj,]
b3$ndh = as.numeric(b3$ndh)
b3$xmt = as.numeric(b3$xmt)

b5 = read.table('data/itc-nacl/B5.txt',sep='\t')
colnames(b5) = c('dh','injv','xt','mt','xmt','ndh','fit','dy')
b5 = b5[2:ninj,]
b5$ndh = as.numeric(b5$ndh)
b5$xmt = as.numeric(b5$xmt)

b7 = read.table('data/itc-nacl/B7.txt',sep='\t')
colnames(b7) = c('dh','injv','xt','mt','xmt','ndh')
b7 = b7[2:ninj,]
b7$ndh = as.numeric(b7$ndh)
b7$xmt = as.numeric(b7$xmt)

nacl_params = data.frame(
  color=c('#fe9929','#ec7014','#cc4c02','#8c2d04'),
  conc=c(0,50,150,500),
  expta=c('a1','a3','a5','a7'),
  exptb=c('b1','b3','b5','b7')
)


plot(NA, NA, xlim=c(0,.50), ylim=c(-1e5, 2e4), axes=F, ann=F, xaxs='i', yaxs='i')
abline(h=0, lwd=0.25)
axis(side=1, at=(0:5)/10)
mtext(side=1, line=2.5, text='molar ratio')
axis(side=2, at=c(-1e5, -5e4, 0, 5e4), labels=c('-100', '-50', '0', ''), las=2)
mtext(side=2, line=2.5, text='kcal/mol injectant')

for (i in c(1,3,5,7)) {
  dataname = paste('a',i,sep='')
  dataset = get(dataname)
  color = nacl_params$color[nacl_params$expta==dataname]
  points(dataset$xmt, dataset$ndh, pch=15, cex=0.6, col=color)
  points(dataset$xmt,dataset$fit,type='l',lwd=1,col=color)
  if (i == 7) {
    segments(x0=min(dataset$xmt),x1=max(dataset$xmt),y0=mean(dataset$ndh), col=color)
  }
}

legend(x=0.35,y=-5.5e4,legend=nacl_params$conc,col=nacl_params$color,pch=15,lwd=1,text.col=nacl_params$color,title.col='black',title='[NaCl] (mM)',bty='n',cex=0.8)
mtext('A', side=3, cex=2, adj = -0.1, line = 0.5)

plot(NA, NA, xlim=c(0,.50), ylim=c(-1e5, 2e4), axes=F, ann=F, xaxs='i', yaxs='i')
abline(h=0, lwd=0.25)
axis(side=1, at=(0:5)/10)
mtext(side=1, line=2.5, text='molar ratio')
axis(side=2, at=c(-1e5, -5e4, 0, 5e4), labels=c('-100', '-50', '0', ''), las=2)
mtext(side=2, line=2.5, text='kcal/mol injectant')

for (i in c(1,3,5,7)) {
  dataname = paste('b',i,sep='')
  dataset = get(dataname)
  color = nacl_params$color[nacl_params$exptb==dataname]
  points(dataset$xmt, dataset$ndh, pch=15, cex=0.6, col=color)
  points(dataset$xmt,dataset$fit,type='l',lwd=1,col=color)
  if (i == 7) {
    segments(x0=min(dataset$xmt),x1=max(dataset$xmt),y0=mean(dataset$ndh), col=color)
  }
}

legend(x=0.35,y=-5.5e4,nacl_params$conc,col=nacl_params$color,pch=15,lwd=1,text.col=nacl_params$color,title.col='black',title='[NaCl] (mM)',bty='n',cex=0.8)
mtext('B', side=3, cex=2, adj = -0.1, line = 0.5)

dev.off() #### --- END FIGURE 3









#### FIGURE 4: ACTIVITY IN SCN2A CELLS
imgsave(paste('figures/figure-4.',imgmode,sep=''),width=6.5*resx,height=3*resx,res=resx)

plot_params = read.table('data/plot_params.tsv',sep='\t',header=T,comment.char='')
dat = read.table('data/scn2a/scn2a_doseresp.tsv',sep='\t',header=T)
dat$color = plot_params$color[match(dat$treatment,plot_params$treatment)]

leg = plot_params[plot_params$treatment %in% dat$treatment,c('treatment','color')]
leg = leg[with(leg, order(treatment)),]

xmin = 0.015
xzero = 0.026
xmax = 25
xlim = c(xmin,xmax)
ylim = c(0, 1.5)

doses = data.frame(dose=sort(unique(dat$dose)))
doses$xplot = log10(doses$dose)
doses$xplot[doses$dose==0] = log10(xzero)
doses$disp = trimws(formatC(doses$dose, format='fg', digits=2))

par(mfrow=c(1,3), mar=c(4, 5, 3, 2))

set.seed(1)
dat$xplot = jitter(log10(dat$dose), .75)
dat$xplot[dat$dose==0] = jitter(log10(rep(xzero, sum(dat$dose==0))), .75)

panel = 1

for (var in unique(dat$measurement)) {
  plot(NA, NA, xlim=log10(xlim), ylim=ylim, xaxs='i', yaxs='i', ann=F, axes=F)
  axis(side=1, at=log10(xlim), labels=NA, lwd=1, lwd.ticks=0)
  axis(side=1, at=doses$xplot, labels=NA, lwd=0, lwd.ticks=1, tck=-0.01)
  par(xpd=T)
  text(x=doses$xplot, y=rep(-.1, nrow(doses)), labels=doses$disp, srt=45)
  par(xpd=F)
  axis.break(axis=1, breakpos=log10(xzero)+.25, style = 'slash', brw=0.05)
  mtext(side=1, line=2.5, text='dose (µM)')
  axis(side=2, at=(0:6)/4, labels=percent((0:6)/4), las=2)
  mtext(side=2, line=3.5, text='relative amount (% saline)')
  mtext(side=3, line=0.5, text=var)
  abline(h=1, lwd=0.5, lty=3)
  subs = subset(dat, measurement==var)
  points(subs$xplot, subs$pct_pbs, col=subs$color, pch=20, cex=0.6)
  for (aso in unique(subs$treatment)) {
    dr = subset(subs, treatment==aso)
    
    m = drm(pct_pbs ~ dose, data=dr, fct=LL.4())
    x = (0:2000)/100
    y = predict(m, newdata=data.frame(dose=x))
    
    points(log10(x), y, type='l', col=plot_params$color[plot_params$treatment==aso])
  }
  mtext(LETTERS[panel], side=3, cex=2, adj = -0.3, line = 0.5)
  if (panel == 1) {
    legend('bottomleft',leg$treatment,col=leg$color,text.col=leg$color,pch=20,lwd=1,bty='n')
  }
  panel = panel + 1
}

dev.off() #### --- END FIGURE 4













#### FIGURE 5: ELISA MEASUREMENT OF CSF PrP
imgsave(paste('figures/figure-5.',imgmode,sep=''),width=3.25*resx,height=5*resx,res=resx)

par(mfrow=c(3,1), mar=c(5,4,3,2))

xlim = c(-10.5, -2.5) # use broken axis to set -10 = no ASO, -9 = 1ng/mL, -3 = 1 mg/mL etc.
ylim = c(0, 400)
xbigs = -10:-3
xbigs_labs1 = c('0','1','10','100','1','10','100','1')
xbigs_labs2 = c('','ng/mL','','','µg/mL','','','mg/mL')
xats = log10(rep(1:9,7) * rep(10^(-9:-3),each=9))

p44 = read.table('data/elisa/44_summary.tsv',sep='\t',header=T)
p51 = read.table('data/elisa/51_summary.tsv',sep='\t',header=T)

p51$x = c(-3, -5, -7, -8, -9, -10, NA)
p51$l95 = p51$ngml_av - 1.96 * p51$se_mean
p51$u95 = p51$ngml_av + 1.96 * p51$se_mean

# read p44 in raw again - use only the 1:50 dilution so everything is as comparable as possible
# and average the two control aliquots together
p44_raw = read.table('data/elisa/44.tsv',sep='\t',header=T)
p44_raw$sample = p44_raw$detail
p44_raw$sample[p44_raw$detail %in% c('9-IPC-control-A-0-ng/mL','9-IPC-control-B-0-ng/mL')] = '9-IPC-control-0-ng/mL'
p44 = sqldf("
            select   sample, avg(ngml_trunc) ngml_av, stdev(ngml_trunc)/sqrt(count(*)) se_mean
            from     p44_raw
            where    stype = 'sample'
            and      detail like '%9-IPC%'
            and      length(detail) > 5
            and      dilution = 50
            group by 1
            order by 1
            ;")
p44$x = c(-3, -9, -8, -5, -7, -10, -3, -9, -8, -5, -7)
p44$l95 = p44$ngml_av - 1.96 * p44$se_mean
p44$u95 = p44$ngml_av + 1.96 * p44$se_mean

aso = subset(p44, grepl('ASO|control', sample))
plot(NA, NA, xlim=xlim, ylim=ylim, axes=F, ann=F, xaxs='i', yaxs='i')
axis(side=2, at=(0:4)*100, las=2)
mtext(side=2, line=2.5, text='[PrP] (ng/mL)')
axis(side=1, at=xlim, labels=NA, lwd=1, lwd.ticks=0)
axis(side=1, at=xats, labels=NA, lwd=0, lwd.ticks=1, tck=-0.025)
axis(side=1, at=xbigs, labels=xbigs_labs1, lwd=0, lwd.ticks=1, tck=-0.05)
axis(side=1, at=xbigs, labels=xbigs_labs2, line=1, lwd=0, lwd.ticks=0, cex.axis=0.9)
axis.break(axis=1,breakpos=-9.5,style = 'slash',brw=0.05)
mtext(side=1, line=3.5, text='[ASO]')
segments(x0=aso$x, y0=aso$l95, y1=aso$u95)
points(x=aso$x, y=aso$ngml_av, pch=19)
mtext('A', side=3, cex=2, adj = -0.1, line = 0.5)

hep = subset(p44, grepl('heparin|control', sample))
plot(NA, NA, xlim=xlim, ylim=ylim, axes=F, ann=F, xaxs='i', yaxs='i')
axis(side=2, at=(0:4)*100, las=2)
mtext(side=2, line=2.5, text='[PrP] (ng/mL)')
axis(side=1, at=xlim, labels=NA, lwd=1, lwd.ticks=0)
axis(side=1, at=xats, labels=NA, lwd=0, lwd.ticks=1, tck=-0.025)
axis(side=1, at=xbigs, labels=xbigs_labs1, lwd=0, lwd.ticks=1, tck=-0.05)
axis(side=1, at=xbigs, labels=xbigs_labs2, line=1, lwd=0, lwd.ticks=0, cex.axis=0.9)
axis.break(axis=1,breakpos=-9.5,style = 'slash',brw=0.05)
mtext(side=1, line=3.5, text='[heparin]')
segments(x0=hep$x, y0=hep$l95, y1=hep$u95)
points(x=hep$x, y=hep$ngml_av, pch=19)
mtext('B', side=3, cex=2, adj = -0.1, line = 0.5)

ylim = c(0,200)

plot(NA, NA, xlim=xlim, ylim=ylim, axes=F, ann=F, xaxs='i', yaxs='i')
axis(side=2, at=(0:4)*50, las=2)
mtext(side=2, line=2.5, text='[PrP] (ng/mL)')
axis(side=1, at=xlim, labels=NA)
axis(side=1, at=xats, labels=NA, lwd=0, lwd.ticks=1, tck=-0.025)
axis(side=1, at=xbigs, labels=xbigs_labs1, lwd=0, lwd.ticks=1, tck=-0.05)
axis(side=1, at=xbigs, labels=xbigs_labs2, line=0.75, lwd=0, lwd.ticks=0)
axis.break(axis=1,breakpos=-9.5,style = 'slash',brw=0.05)
mtext(side=1, line=3.5, text='[ASO]')
segments(x0=p51$x, y0=p51$l95, y1=p51$u95)
points(x=p51$x, y=p51$ngml_av, pch=19)
mtext('C', side=3, cex=2, adj = -0.1, line = 0.5)

dev.off()




