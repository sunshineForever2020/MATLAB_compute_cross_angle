function anglertm = compute_angle_rtm(tazi, telev, trange)
% 该函数计算天线指向向量与探测点处磁力线向量的夹角
% 问题归约为sysir相位中心点r，探测点t，t点处磁场强度点m三者构成的角rtm
% 输入： 
%     tazi    : 探测点在雷达坐标系内的方位角  单位 度   标量 | 向量
%     telev   : 探测点在雷达坐标系内的俯仰角  单位 度   标量 | 向量
%     trange  : 探测点在雷达坐标系内的斜距离  单位 度   标量 | 向量
% 输出：
%     anglertm: 天线指向rt与t点磁场强度tm的夹角 单位度 标量 | 向量
narginchk(3,3);
% 地球模型参考椭球， wgs84, 单位米
wgs84 = wgs84Ellipsoid('meter');
% sysir 雷达相位中心地理坐标（已知）
rlat0 = dms2degrees([18 20 56.79]);
rlon0 = dms2degrees([109 37 19.86]);
rh0 = 54;
% 计算雷达相位中心地心地固坐标
[rx0, ry0, rz0] = geodetic2ecef(wgs84, rlat0, rlon0, rh0);
rxyz0 = [rx0 ry0 rz0]; 
% 计算目标的地理坐标
[tlat, tlon, th] = aer2geodetic(tazi, telev, trange, rlat0, rlon0, rh0, wgs84);
% 计算相位中心到目标的指向，即天线指向向量rt
[tx, ty, tz] = geodetic2ecef(wgs84, tlat, tlon, th);
txyz = [tx ty tz]; 
rt = txyz - rxyz0; 
% 计算目标到磁场的向量tm 
time = repmat(decyear(now), length(tazi), 1);
mned = igrfmagm(th, tlat, tlon, time, 13);
mnorth = mned(:,1); meast = mned(:, 2); mdown = mned(:, 3);
[mx, my, mz] = ned2ecef(mnorth, meast, mdown, tlat, tlon, th, wgs84);
mxyz = [mx my mz]; 
tm = mxyz - txyz;
% 计算角rtm，即天线指向与目标点处磁力线夹角
anglertm = rad2deg(acos(dot(rt, tm, 2) ./ ...
    (vecnorm(rt, 2, 2) .* vecnorm(tm, 2, 2))));
end

