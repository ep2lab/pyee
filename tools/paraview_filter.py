data = inputs[0]

Ezr = data.PointData["real_Ez"]
Ezi = data.PointData["imag_Ez"]
output.PointData.append((Ezi[:]**2 + Ezr[:]**2)**0.5,"Ez_mag")
output.PointData.append(arctan2(Ezi[:],Ezr[:]),"Ez_phs")

Exr = data.PointData["real_Ex"]
Exi = data.PointData["imag_Ex"]
output.PointData.append((Exi[:]**2 + Exr[:]**2)**0.5,"Ex_mag")
output.PointData.append(arctan2(Exi[:],Exr[:]),"Ex_phs")

Eyr = data.PointData["real_Ey"]
Eyi = data.PointData["imag_Ey"]
output.PointData.append((Eyi[:]**2 + Eyr[:]**2)**0.5,"Ey_mag")
output.PointData.append(arctan2(Eyi[:],Eyr[:]),"Ey_phs")