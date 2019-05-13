load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

begin

; Parameters:
; -----------

datadir = getenv("datadir")
filebase = getenv("filebase")
date = getenv("date")
difference = getenv("difference")
if (difference.eq."true")
  filebase_ref = getenv("filebase_ref")
  date_ref = getenv("date_ref")
end if

; Load data:
; ----------

data = addfile(datadir + "/" + filebase + "." + date + ".nc","r")
lx = isfilevar(data,"x")
lq = isfilevar(data,"q")
lw = isfilevar(data,"u").and.isfilevar(data,"v")
if (lx) then
   x = data->x
   dims = dimsizes(x)
end if
if (lq) then
   q = data->q
   dims = dimsizes(q)
end if
if (lw) then
   u = data->u
   v = data->v
   dims = dimsizes(u)
end if
if (difference.eq."true")
   data_ref = addfile(datadir + "/" + filebase_ref + "." + date_ref + ".nc","r")
   if (lx) then
      x = x-data_ref->x
   end if
   if (lq) then
      q = q-data_ref->q
   end if
   if (lw) then
      u = u-data_ref->u
      v = v-data_ref->v
   end if
end if

; Resources:
; ----------

res = True
res@gsnDraw = False 
res@gsnFrame = False 
res@gsnMaximize = True
res@gsnAddCyclic = False

res@cnFillOn = True
res@cnFillMode = "AreaFill"
res@trGridType = "TriangularMesh"
res@cnMonoFillPattern = True
res@cnMonoFillColor = False
res@lbLabelBarOn = False
res@cnInfoLabelOn = False
res@cnLineLabelsOn = False
res@cnLinesOn = False
res@cnNoDataLabelOn = False
res@cnMissingValFillColor = -1
res@vpHeightF = 0.1
res@vpWidthF = 0.2
res@tmXBOn = False
res@tmXTOn = False
res@tmYROn = False
res@tmYLOn = False

vres = True
vres@gsnDraw = False 
vres@gsnFrame = False 
vres@gsnMaximize = True
vres@gsnAddCyclic = False
vres@vcGlyphStyle = "CurlyVector"
vres@vcLineArrowThicknessF = 1.0
vres@vcRefLengthF = 0.04
vres@vcRefAnnoOn = False
vres@lbLabelBarOn = False
vres@lbOrientation = "vertical"
vres@vcMinDistanceF = 0.01
vres@vcRefMagnitudeF = 10.0
vres@vcLineArrowHeadMaxSizeF = 0.003
vres@vcLineArrowThicknessF = 2.0
vres@vpHeightF = 0.1
vres@vpWidthF = 0.2
vres@tmXBOn = False
vres@tmXTOn = False
vres@tmYROn = False
vres@tmYLOn = False

txres = True
txres@txFontHeightF = 0.012
txres@txJust = "CenterCenter"

resP = True
resP@gsnFrame = False
resP@gsnPanelYWhiteSpacePercent = 3
resP@gsnPanelXWhiteSpacePercent = 3
resP@gsnPanelTop = 0.9
resP@gsnPanelLeft = 0.03

; Display:
; --------

system("mkdir -p " + datadir + "/fig")

nplot = 0
if (lx) then
   nplot = nplot+1
end if
if (lq) then
   nplot = nplot+1
end if
if (lw) then
   nplot = nplot+1
end if
plot = new(nplot*dims(0),graphic)
wks_type = "png"
wks_type@wkWidth = 2500
wks_type@wkHeight = 2500 

if (difference.eq."true")
   output = datadir + "/fig/" + filebase + "-" + filebase_ref + "." + date
else
   output = datadir + "/fig/" + filebase + "." + date
end if
wks = gsn_open_wks(wks_type,output)
if (difference.eq."true") then
   gsn_define_colormap(wks,"BlueWhiteOrangeRed")
   res@cnLevelSelectionMode = "ManualLevels"
else
   gsn_define_colormap(wks,"WhiteBlueGreenYellowRed")
end if
iplot = 0
do iz=0,dims(0)-1
   if (lx) then
      if (difference.eq."true") then
         bound = max(abs(x(iz,:,:)))
         res@cnMaxLevelValF = bound
         res@cnMinLevelValF = -bound
         res@cnLevelSpacingF = 0.05*bound
      end if
      res@tiMainString = "Streamfunction @ level " + iz
      plot(iplot) = gsn_csm_contour(wks,x(iz,:,:),res)
      iplot = iplot+1
   end if
   if (lq) then
      if (difference.eq."true") then
         bound = max(abs(q(iz,:,:)))
         res@cnMaxLevelValF = bound
         res@cnMinLevelValF = -bound
         res@cnLevelSpacingF = 0.05*bound
      end if
      res@tiMainString = "Potential vorticity @ level " + iz
      plot(iplot) = gsn_csm_contour(wks,q(iz,:,:),res)
      iplot = iplot+1
   end if
   if (lw) then
      vres@tiMainString = "Wind - level " + iz
      if (difference.eq."true") then
         plot(iplot) = gsn_csm_vector(wks,u(iz,:,:),v(iz,:,:),vres)
      else
         plot(iplot) = gsn_csm_vector_scalar(wks,u(iz,:,:),v(iz,:,:),sqrt(u(iz,:,:)^2+v(iz,:,:)^2),vres)
      end if
      iplot = iplot+1
   end if
end do

gsn_panel(wks,plot,(/dims(0),nplot/),resP) 
frame(wks)
delete(wks)

system("mogrify -trim " + output + ".png")

end
