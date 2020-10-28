load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

; Parameters:
; -----------


datadir = getenv("datadir")
filebase = getenv("filebase")

; Load data and get sizes:
; ------------------------

Jb = readAsciiTable("tmp/Jb",1,"float",0)
Jo_vec = readAsciiTable("tmp/Jo",1,"float",0)
Jo_n_vec = readAsciiTable("tmp/Jo_n",1,"float",0)
dims = dimsizes(Jb)
nouter = dims(0)-1

obstypes = readAsciiTable("tmp/obstypes",1,"string",0)
dims = dimsizes(obstypes)
nobstypes = dims(0)

if (numAsciiRow("tmp/Jc").gt.0) then
   Jc_vec = readAsciiTable("tmp/Jc",1,"float",0)
   dims = dimsizes(Jc_vec)
   nts = dims(0)/(nouter+1)
else
   Jc_vec = new((/nouter+1,1/),"float")
   Jc_vec@_FillValue = Jb@_FillValue
   Jc_vec = Jc_vec@_FillValue
   nts = 1
end if

cost_quad = (numAsciiRow("tmp/Jb_quad").gt.0)
if (cost_quad) then
   ninner = readAsciiTable("tmp/ninner",1,"integer",0)
   Jb_quad = readAsciiTable("tmp/Jb_quad",1,"float",0)
   JoJc_quad = readAsciiTable("tmp/JoJc_quad",1,"float",0)
   dims = dimsizes(Jb_quad)
   ninner_max = max(ninner)
end if

; Prepare data:
; -------------

J = new((/4,nouter+1/),float)
Jo_n = new((/nobstypes,nouter+1/),float)
iotot = 0
itstot = 0
do iouter=0,nouter
   J(1,iouter) = 0.0
   do its=0,nts-1
      J(1,iouter) = J(1,iouter)+Jc_vec(itstot,0)
      itstot = itstot+1
   end do
   J(2,iouter) = 0.0
   do iobstypes=0,nobstypes-1
      J(2,iouter) = J(2,iouter)+Jo_vec(iotot,0)
      Jo_n(nobstypes-1-iobstypes,iouter) = Jo_n_vec(iotot,0)
      iotot = iotot+1
   end do
   J(3,iouter) = Jb(iouter,0)
   J(0,iouter) = sum(J(1:3,iouter))
end do
obstypes(:,0) = obstypes(::-1,0)
xouter = fspan(0,nouter,nouter+1)

if (cost_quad) then
   J_quad = new((/3,ninner_max,nouter/),float)
   itot = 0
   do iouter=0,nouter-1
      do iinner=0,ninner(iouter,0)-1
         J_quad(1,iinner,iouter) = JoJc_quad(itot,0)
         J_quad(2,iinner,iouter) = Jb_quad(itot,0)
         J_quad(0,iinner,iouter) = sum(J_quad(1:2,iinner,iouter))
         itot = itot+1
      end do
   end do
   xinner = fspan(1,ninner_max,ninner_max)
end if

; Resources :
; -----------

xyres = True
xyres@gsnFrame = False
xyres@gsnDraw = False
xyres@tmXBLabelFontHeightF = 0.015
xyres@tmYLLabelFontHeightF = 0.015
xyres@xyLineThicknessF = 10.0
xyres@xyCurveDrawOrder = "PreDraw"
xyres@xyDashPattern = 0
xyres@xyMarkLineMode = "MarkLines"
xyres@xyMarker = 16
xyres@pmLegendDisplayMode = "Always"
xyres@tiMainFontHeightF = 0.02
xyres@tiXAxisFontHeightF = 0.02
xyres@tiYAxisFontHeightF = 0.02
xyres@vpHeightF = 0.3
xyres@vpWidthF = 0.4
xyres@lgJustification = "TopRight"
xyres@pmLegendSide = "Top"
xyres@pmLegendOrthogonalPosF = -0.503
xyres@pmLegendParallelPosF = 1.0
xyres@pmLegendWidthF = 0.07
xyres@pmLegendHeightF = 0.1
xyres@lgBoxMinorExtentF = 0.4
xyres@lgLabelFontHeightF = 0.02

txres = True
txres@txFontHeightF = 0.02
txres@txJust = "CenterCenter"

resP = True
resP@gsnFrame = False
resP@gsnPanelYWhiteSpacePercent = 5
resP@gsnPanelXWhiteSpacePercent = 5
;resP@gsnPanelTop = 0.9
;resP@gsnPanelLeft = 0.03

; Display :
; ---------

system("mkdir -p " + datadir + "/fig")

plot = new(2,graphic)

output = datadir + "/fig/" + filebase + "_cost"
wks_type = "png"
wks_type@wkWidth = 2500
wks_type@wkHeight = 2500
wks = gsn_open_wks(wks_type,output)
gsn_define_colormap(wks,"BlAqGrYeOrReVi200")

xyres@tiMainString = "Nonlinear cost function"
xyres@tiXAxisString = "Outer iteration"
xyres@trXMinF = 0
xyres@trXMaxF = nouter
xyres@tmXBMinorOn = False
xyres@tmXTMinorOn = False
xyres@tmXBMode = "Explicit"
xyres@tmXBValues = toint(xouter)
xyres@tmXBLabels = xyres@tmXBValues
xyres@xyYStyle = "Linear"
xyres@trYMinF = 0.0
xyres@trYMaxF = 1.1*max(J)
xyres@xyLineColors = (/1,140,2,170/)
xyres@xyMarkerColors = xyres@xyLineColors
xyres@xyExplicitLegendLabels = "  " + (/"J ","J~B~c~N~","J~B~o~N~","J~B~b~N~"/)
plot(0) = gsn_csm_xy(wks,xouter,J,xyres)

delete(xyres@xyLineColors)
delete(xyres@xyMarkerColors)
delete(xyres@xyExplicitLegendLabels)

xyres@tiMainString = "Nonlinear J~B~o~N~ / n"
xyres@xyYStyle = "Log"
xyres@trYMinF = 0.5*min(Jo_n)
xyres@trYMaxF = 2.0*max(Jo_n)
xyres@xyLineColors = toint(fspan(2,90,nobstypes))
xyres@xyMarkerColors = xyres@xyLineColors
xyres@xyExplicitLegendLabels = "  " + obstypes(:,0)
plot(1) = gsn_csm_xy(wks,xouter,Jo_n,xyres)

gsn_text_ndc(wks,filebase + " with " + nouter + " outer iterations",0.52,0.75,txres)

gsn_panel(wks,plot,(/1,2/),resP)
frame(wks)
delete(wks)

system("mogrify -trim " + output + ".png")

if (cost_quad) then
   delete(xyres@xyLineColors)
   delete(xyres@xyMarkerColors)
   delete(xyres@tmXBValues)
   delete(xyres@tmXBLabels)
   delete(plot)

   plot = new(nouter,graphic)

   output = datadir + "/fig/" + filebase + "_cost_quad"
   wks_type = "png"
   wks_type@wkWidth = 2500
   wks_type@wkHeight = 2500
   wks = gsn_open_wks(wks_type,output)
   gsn_define_colormap(wks,"BlAqGrYeOrReVi200")

   xyres@tiXAxisString = "Inner iteration"
   xyres@trXMinF = 1
   xyres@tmXBMinorOn = False
   xyres@tmXTMinorOn = False
   xyres@tmXBMode = "Explicit"
   xyres@tmXBValues = toint(xinner)
   xyres@xyYStyle = "Linear"
   xyres@xyLineColors = (/1,2,170/)
   xyres@xyMarkerColors = xyres@xyLineColors
   xyres@xyExplicitLegendLabels = "  " + (/"J ","J~B~o~N~+J~B~c~N~","J~B~b~N~"/)

   do iouter=0,nouter-1
      xyres@tiMainString = "Quad. cost function, outer iteration " + (iouter+1)
      xyres@trXMaxF = ninner(iouter,0)

      xyres@tmXBLabels = xyres@tmXBValues
      if (min(J_quad(:,0:ninner(iouter,0)-1,iouter)).gt.0) then
         xyres@trYMinF = 0.9*min(J_quad(:,0:ninner(iouter,0)-1,iouter))
      else
         xyres@trYMinF = 1.1*min(J_quad(:,0:ninner(iouter,0)-1,iouter))
      end if
      xyres@trYMaxF = 1.1*max(J_quad(:,0:ninner(iouter,0)-1,iouter))
      plot(iouter) = gsn_csm_xy(wks,xinner(0:ninner(iouter,0)-1),J_quad(:,0:ninner(iouter,0)-1,iouter),xyres)
   end do

   gsn_text_ndc(wks,filebase + " with " + nouter + " outer iterations",0.52,0.75,txres)

   gsn_panel(wks,plot,(/1,nouter/),resP)
   frame(wks)
   delete(wks)

   system("mogrify -trim " + output + ".png")
end if

end
