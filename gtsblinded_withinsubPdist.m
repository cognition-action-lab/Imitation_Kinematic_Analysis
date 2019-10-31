%quick within-subject analysis of gesture-to-sight data
%for all selected files, go through and run a procrustes comparison of
%blinded vs non-blinded
%compile the pdist values

clear all;

[fnames,fpath] = uigetfile('*markdata.mat','Select data file to process.','Multiselect','on');

if ~iscell(fnames)
    fnames = {fnames};
end

xlsfiles={fnames};
[~,idx]=sort(fnames);
fnames=fnames(idx);
clear xlsfiles;

%function to calculate the rotation matrix given the Euler angles
REuler = @(e,a,r) [  cosd(e)*cosd(a)                              cosd(e)*sind(a)                              -sind(e);
                   -(cosd(r)*sind(a))+(sind(r)*sind(e)*cosd(a))   (cosd(r)*cosd(a))+(sind(r)*sind(e)*sind(a))   sind(r)*cosd(e);
                    (sind(r)*sind(a))+(cosd(r)*sind(e)*cosd(a))  -(sind(r)*cosd(a))+(cosd(r)*sind(e)*sind(a))   cosd(r)*cosd(e)];



blk = 10;

for a = 1:length(fnames)
    
    tmp = load(fullfile(fpath,fnames{a}));
    
    for c = 1:41 %size(tmp.BlockData{blk}.pos,1)
        %pull the position data
        X = tmp.BlockData{blk}.pos{c,1};
        Y = tmp.BlockData{blk}.pos{c,2};
        
        if isempty(X) || isempty(Y) || any(isnan(tmp.BlockData{blk}.inds(c,1).Full)) || any(isnan(tmp.BlockData{blk}.inds(c,2).Full))
            PDist.Pos(a,c) = NaN;
            PDist.ArmPlane(a,c) = NaN;
            PDist.JA(a,c) = NaN;
            PDist.Vel(a,c) = NaN;
            PDist.Vel3D(a,c) = NaN;
            continue;
        end
        
        %extract out the good part of the reach to analyze
        X = X(tmp.BlockData{blk}.inds(c,1).Full(1):tmp.BlockData{blk}.inds(c,1).Full(2),:,:);
        for d = 1:size(X,2)
            for e = 1:size(X,3)
                X(:,d,e) = sgolayfilt(X(:,d,e),2,min([101,2*floor((size(X,1)-2)/2)+1])); %smooth the data channels
            end
        end
        Y = Y(tmp.BlockData{blk}.inds(c,2).Full(1):tmp.BlockData{blk}.inds(c,2).Full(2),:,:);
        for d = 1:size(Y,2)
            for e = 1:size(Y,3)
                Y(:,d,e) = sgolayfilt(Y(:,d,e),2,min([101,2*floor((size(Y,1)-2)/2)+1])); %smooth the data channels
            end
        end
        
        %pull the rotation matrix
        Xrotang = tmp.BlockData{blk}.rotang{c,1};
        Yrotang = tmp.BlockData{blk}.rotang{c,2};
        Xrotang = Xrotang(:,:,tmp.BlockData{blk}.inds(c,1).Full(1):tmp.BlockData{blk}.inds(c,1).Full(2),:);
        for d = 1:3
            for e = 1:3
                for f = 1:size(Xrotang,4)
                    Xrotang(d,e,:,f) = sgolayfilt(Xrotang(d,e,:,f),2,min([101,2*floor((size(Xrotang,3)-2)/2)+1])); %smooth the data channels
                end
            end
        end
        Yrotang = Yrotang(:,:,tmp.BlockData{blk}.inds(c,2).Full(1):tmp.BlockData{blk}.inds(c,2).Full(2),:);
        for d = 1:3
            for e = 1:3
                for f = 1:size(Yrotang,4)
                    Yrotang(d,e,:,f) = sgolayfilt(Yrotang(d,e,:,f),2,min([101,2*floor((size(Yrotang,3)-2)/2)+1])); %smooth the data channels
                end
            end
        end
        
        %pull the joint-angle data
        Xja = tmp.BlockData{blk}.ja{c,1}';
        Yja = tmp.BlockData{blk}.ja{c,2}';
        Xja = Xja(tmp.BlockData{blk}.inds(c,1).Full(1):tmp.BlockData{blk}.inds(c,1).Full(2),:);
        for d = 1:size(Xja,2)
            Xja(:,d) = sgolayfilt(Xja(:,d),2,min([101,2*floor((size(Xja,1)-2)/2)+1])); %smooth the data channels
        end
        Yja = Yja(tmp.BlockData{blk}.inds(c,2).Full(1):tmp.BlockData{blk}.inds(c,2).Full(2),:);
        for d = 1:size(Xja,2)
            Yja(:,d) = sgolayfilt(Yja(:,d),2,min([101,2*floor((size(Yja,1)-2)/2)+1])); %smooth the data channels
        end
        
        %pull out the velocity data
        Xv = tmp.BlockData{blk}.vel{c,1};
        Yv = tmp.BlockData{blk}.vel{c,2};
        Xv = Xv(tmp.BlockData{blk}.inds(c,1).Full(1):tmp.BlockData{blk}.inds(c,1).Full(2),:,:);
        for d = 1:size(Xv,2)
            for e = 1:size(Xv,3)
                Xv(:,d,e) = sgolayfilt(Xv(:,d,e),2,min([101,2*floor((size(Xv,1)-2)/2)+1])); %smooth the data channels
            end
        end
        Yv = Yv(tmp.BlockData{blk}.inds(c,2).Full(1):tmp.BlockData{blk}.inds(c,2).Full(2),:,:);
        for d = 1:size(Xv,2)
            for e = 1:size(Xv,3)
                Yv(:,d,e) = sgolayfilt(Yv(:,d,e),2,min([101,2*floor((size(Yv,1)-2)/2)+1])); %smooth the data channels
            end
        end
        Xv3d = squeeze(sqrt(sum(Xv.^2,2)));
        Yv3d = squeeze(sqrt(sum(Yv.^2,2)));
        
        
        %time-normalize the trajectories
        N = 1000;
        
        t = tmp.BlockData{blk}.time{c,1}(tmp.BlockData{blk}.inds(c,1).Full(1):tmp.BlockData{blk}.inds(c,1).Full(2))';
        t = t-t(1);
        dur = t(end);
        dN = dur/N;
        tinterp = [0:dN:dur];
        
        Xrsmp = permute(csapi(t,permute(X,[2,3,1]),tinterp),[3,1,2]);  %each column of X is a sample
        Xrsmpreshape = [Xrsmp(:,:,1) Xrsmp(:,:,2) Xrsmp(:,:,3) Xrsmp(:,:,4) Xrsmp(:,:,5) Xrsmp(:,:,6) Xrsmp(:,:,7) Xrsmp(:,:,8)];
        Xrotangrsmp = permute(csapi(t,permute(Xrotang,[1,2,4,3]),tinterp),[1,2,4,3]);  %each column of X is a sample
        Xjarsmp = csapi(t,Xja',tinterp)';
        Xvrsmp = permute(csapi(t,permute(Xv,[2,3,1]),tinterp),[3 1 2]);
        Xvrsmpreshape = [Xvrsmp(:,:,1) Xvrsmp(:,:,2) Xvrsmp(:,:,3) Xvrsmp(:,:,4) Xvrsmp(:,:,5) Xvrsmp(:,:,6) Xvrsmp(:,:,7) Xvrsmp(:,:,8)];
        Xv3drsmp = csapi(t,Xv3d',tinterp)';
        
        
        t = tmp.BlockData{blk}.time{c,2}(tmp.BlockData{blk}.inds(c,2).Full(1):tmp.BlockData{blk}.inds(c,2).Full(2))';
        t = t-t(1);
        dur = t(end);
        dN = dur/N;
        tinterp = [0:dN:dur];
        
        Yrsmp = permute(csapi(t,permute(Y,[2,3,1]),tinterp),[3,1,2]);  %each column of X is a sample
        Yrsmpreshape = [Yrsmp(:,:,1) Yrsmp(:,:,2) Yrsmp(:,:,3) Yrsmp(:,:,4) Yrsmp(:,:,5) Yrsmp(:,:,6) Yrsmp(:,:,7) Yrsmp(:,:,8)];
        Yrotangrsmp = permute(csapi(t,permute(Yrotang,[1,2,4,3]),tinterp),[1,2,4,3]);  %each column of X is a sample
        Yjarsmp = csapi(t,Yja',tinterp)';
        Yvrsmp = permute(csapi(t,permute(Yv,[2,3,1]),tinterp),[3 1 2]);
        Yvrsmpreshape = [Yvrsmp(:,:,1) Yvrsmp(:,:,2) Yvrsmp(:,:,3) Yvrsmp(:,:,4) Yvrsmp(:,:,5) Yvrsmp(:,:,6) Yvrsmp(:,:,7) Yvrsmp(:,:,8)];
        Yv3drsmp = csapi(t,Yv3d',tinterp)';
        
        
        %Do a direct Procrustes comparison
        PDist.Pos(a,c) = procrustesMultiScale(Xrsmpreshape,Yrsmpreshape,'reflect',0); %each row is a sample
        
%         figure(101)
%         subplot(6,7,c)
%         Xpca = Xrsmpreshape - mean(Xrsmpreshape);
%         [PCA_coef,PCA_score] = pca(Xpca,'centered',false);
%         Ypca = Yrsmpreshape - mean(Yrsmpreshape);
%         score = Ypca*PCA_coef;
%         
%         plot3(PCA_score(:,1),PCA_score(:,2),PCA_score(:,3),'k-');
%         hold on;
%         plot3(score(:,1),score(:,2),score(:,3),'b-');
%         hold off;
        
        
        %***Compute the arm plane***
        
        %process the X data
        
        %compute the orientation of the shoulder-elbow-wrist plane
        %note: markers are thumb=1, index=2, middle=3, hand=4,
        %                  wrist=5, elbow=6, shoulder=7,
        %                  refshoulder = 8
        vec_shoulder_elbow = Xrsmp(:,:,7) - Xrsmp(:,:,6);
        vec_wrist_elbow = Xrsmp(:,:,5) - Xrsmp(:,:,6);
        n_armplane = cross(vec_shoulder_elbow,vec_wrist_elbow,2);  %this is the z axis
        z_hat = n_armplane./sqrt(sum(n_armplane.^2,2));
        
        %the rotated y axis is along the shoulder-elbow vector at the
        %shoulder
        y_hat = Xrsmp(:,:,6) - Xrsmp(:,:,7);
        y_hat = y_hat./sqrt(sum(y_hat.^2,2));
        
        %the rotated x axis is the cross product of those two vectors
        x_hat = cross(y_hat,z_hat,2);
        
        %compute the rotation matrix for the new coordinate system
        sew_axes = cat(3,x_hat,y_hat,z_hat);  %dim1 = [t]; dim2 = [x0,y0,z0]; dim3 = [x_hat, y_hat, z_hat]
        %this is the transpose of the rotation matrix that takes
        %  standard axes into this new coordinate system
        %  (i.e., IR' = sew_axes...) so we can just back out the three
        %  rotation angles from tmp' (which is tmp for which dims 1 and
        %  2 are swapped)
        temp = permute(sew_axes,[1,3,2]);
        
        ang_armplane_theta = atan2d(temp(:,2,1),temp(:,1,1));  %rotation about z (~left/right)
        ang_armplane_phi = -asind(temp(:,3,1));               %rotation about y
        ang_armplane_psi = atan2d(temp(:,3,2),temp(:,3,3));    %rotation about z (~up/down)
        
        %we will step through and throw out any discontinuities due to
        %angle changes
        ind = find(abs(diff(ang_armplane_phi)) >= 360-30);
        for i = 1:length(ind)
            sign2pi = sign(ang_armplane_phi(ind(i)+1)-ang_armplane_phi(ind(i)));
            ang_armplane_phi(ind(i):end) = ang_armplane_phi(ind(i):end) - sign2pi*(360);
        end
        if min(ang_armplane_phi) > 180 && max(ang_armplane_phi) <= 360 %recenter the angle around 0
            ang_armplane_phi = ang_armplane_phi-360;
        elseif min(ang_armplane_phi) < -360 && max(ang_armplane_phi) < -180
            ang_armplane_phi = ang_armplane_phi+360;
        end
        ind = find(abs(diff(ang_armplane_theta)) >= 360-30);
        for i = 1:length(ind)
            sign2pi = sign(ang_armplane_theta(ind(i)+1)-ang_armplane_theta(ind(i)));
            ang_armplane_theta(ind(i)+1:end) = ang_armplane_theta(ind(i)+1:end) - sign2pi*(360);
        end
        if min(ang_armplane_theta) > 180 && max(ang_armplane_theta) <= 360 %recenter the angle around 0
            ang_armplane_theta = ang_armplane_theta-360;
        elseif min(ang_armplane_theta) < -360 && max(ang_armplane_theta) < -180
            ang_armplane_theta = ang_armplane_theta+360;
        end
        ind = find(abs(diff(ang_armplane_psi)) >= 360-30);
        for i = 1:length(ind)
            sign2pi = sign(ang_armplane_psi(ind(i)+1)-ang_armplane_psi(ind(i)));
            ang_armplane_psi(ind(i)+1:end) = ang_armplane_psi(ind(i)+1:end) - sign2pi*(360);
        end
        if min(ang_armplane_psi) > 180 && max(ang_armplane_psi) <= 360 %recenter the angle around 0
            ang_armplane_psi = ang_armplane_psi-360;
        elseif min(ang_armplane_psi) < -360 && max(ang_armplane_psi) < -180
            ang_armplane_psi = ang_armplane_psi+360;
        end
        
        
        %elbow angle is the angle betwen the two vectors
        elbow_ang = acosd( dot(vec_shoulder_elbow,vec_wrist_elbow,2) ./ (vecnorm(vec_shoulder_elbow,2,2).*vecnorm(vec_wrist_elbow,2,2)) );
        %throw out any discontinuities in angle changes
        ind = find(abs(diff(elbow_ang)) >= 360-30);
        for i = 1:length(ind)
            sign2pi = sign(elbow_ang(ind(i)+1)-elbow_ang(ind(i)));
            elbow_ang(ind(i):end) = elbow_ang(ind(i):end) - sign2pi*(360);
        end
        if min(elbow_ang) > 180 && max(elbow_ang) <= 360 %recenter the angle around 0
            elbow_ang = elbow_ang-360;
        elseif min(elbow_ang) < -360 && max(elbow_ang) < -180
            elbow_ang = elbow_ang+360;
        end
        
        
        %wrist rotation is the angle between the z axis of the wrist
        %tracker (in the standard coordinte frame) and the normal of the plane
        %for d = 1:size(BlockData{blk}.Grp(grp).mvmt(b,c).limbang.wrist,1)
        for d = 1:size(Xrsmp(:,:,6),1)
            %RotWrist = REuler(BlockData{blk}.Grp(grp).mvmt(b,c).limbang.wrist(d,1),BlockData{blk}.Grp(grp).mvmt(b,c).limbang.wrist(d,2),BlockData{blk}.Grp(grp).mvmt(b,c).limbang.wrist(d,3));
            %z_wrist(d,:) = RotWrist*[0 0 -1]';
            z_wrist(d,:) = Xrotangrsmp(:,:,d,5)*[0 0 1]';  %wrist rotation matrix times z vector at time d
        end
        wrist_ang = acosd( dot(n_armplane,z_wrist,2) ./ (vecnorm(n_armplane,2,2).*vecnorm(z_wrist,2,2)) );
        %throw out any discontinuities in angle changes
        ind = find(abs(diff(wrist_ang)) >= 360-30);
        for i = 1:length(ind)
            sign2pi = sign(wrist_ang(ind(i)+1)-wrist_ang(ind(i)));
            wrist_ang(ind(i):end) = wrist_ang(ind(i):end) - sign2pi*(360);
        end
        if min(wrist_ang) > 180 && max(wrist_ang) <= 360 %recenter the angle around 0
            wrist_ang = wrist_ang-360;
        elseif min(wrist_ang) < -360 && max(wrist_ang) < -180
            wrist_ang = wrist_ang+360;
        end
        
        %for some reason, the wrist angle looks really unstable. Something
        %isn't quite right here...
        
        Xarm = [ang_armplane_phi,ang_armplane_theta,ang_armplane_psi,elbow_ang,wrist_ang];
        
        
        %process the Y data
        %compute the orientation of the shoulder-elbow-wrist plane
        %note: markers are thumb=1, index=2, middle=3, hand=4,
        %                  wrist=5, elbow=6, shoulder=7,
        %                  refshoulder = 8
        vec_shoulder_elbow = Yrsmp(:,:,7) - Yrsmp(:,:,6);
        vec_wrist_elbow = Yrsmp(:,:,5) - Yrsmp(:,:,6);
        n_armplane = cross(vec_shoulder_elbow,vec_wrist_elbow,2);  %this is the z axis
        z_hat = n_armplane./sqrt(sum(n_armplane.^2,2));
        
        %the rotated y axis is along the shoulder-elbow vector at the
        %shoulder
        y_hat = Yrsmp(:,:,6) - Yrsmp(:,:,7);
        y_hat = y_hat./sqrt(sum(y_hat.^2,2));
        
        %the rotated x axis is the cross product of those two vectors
        x_hat = cross(y_hat,z_hat,2);
        
        %compute the rotation matrix for the new coordinate system
        sew_axes = cat(3,x_hat,y_hat,z_hat);  %dim1 = [t]; dim2 = [x0,y0,z0]; dim3 = [x_hat, y_hat, z_hat]
        %this is the transpose of the rotation matrix that takes
        %  standard axes into this new coordinate system
        %  (i.e., IR' = sew_axes...) so we can just back out the three
        %  rotation angles from tmp' (which is tmp for which dims 1 and
        %  2 are swapped)
        temp = permute(sew_axes,[1,3,2]);
        
        ang_armplane_theta = atan2d(temp(:,2,1),temp(:,1,1));  %rotation about z (~left/right)
        ang_armplane_phi = -asind(temp(:,3,1));               %rotation about y
        ang_armplane_psi = atan2d(temp(:,3,2),temp(:,3,3));    %rotation about z (~up/down)
        
        %we will step through and throw out any discontinuities due to
        %angle changes
        ind = find(abs(diff(ang_armplane_phi)) >= 360-30);
        for i = 1:length(ind)
            sign2pi = sign(ang_armplane_phi(ind(i)+1)-ang_armplane_phi(ind(i)));
            ang_armplane_phi(ind(i):end) = ang_armplane_phi(ind(i):end) - sign2pi*(360);
        end
        if min(ang_armplane_phi) > 180 && max(ang_armplane_phi) <= 360 %recenter the angle around 0
            ang_armplane_phi = ang_armplane_phi-360;
        elseif min(ang_armplane_phi) < -360 && max(ang_armplane_phi) < -180
            ang_armplane_phi = ang_armplane_phi+360;
        end
        ind = find(abs(diff(ang_armplane_theta)) >= 360-30);
        for i = 1:length(ind)
            sign2pi = sign(ang_armplane_theta(ind(i)+1)-ang_armplane_theta(ind(i)));
            ang_armplane_theta(ind(i)+1:end) = ang_armplane_theta(ind(i)+1:end) - sign2pi*(360);
        end
        if min(ang_armplane_theta) > 180 && max(ang_armplane_theta) <= 360 %recenter the angle around 0
            ang_armplane_theta = ang_armplane_theta-360;
        elseif min(ang_armplane_theta) < -360 && max(ang_armplane_theta) < -180
            ang_armplane_theta = ang_armplane_theta+360;
        end
        ind = find(abs(diff(ang_armplane_psi)) >= 360-30);
        for i = 1:length(ind)
            sign2pi = sign(ang_armplane_psi(ind(i)+1)-ang_armplane_psi(ind(i)));
            ang_armplane_psi(ind(i)+1:end) = ang_armplane_psi(ind(i)+1:end) - sign2pi*(360);
        end
        if min(ang_armplane_psi) > 180 && max(ang_armplane_psi) <= 360 %recenter the angle around 0
            ang_armplane_psi = ang_armplane_psi-360;
        elseif min(ang_armplane_psi) < -360 && max(ang_armplane_psi) < -180
            ang_armplane_psi = ang_armplane_psi+360;
        end
        
        
        %elbow angle is the angle betwen the two vectors
        elbow_ang = acosd( dot(vec_shoulder_elbow,vec_wrist_elbow,2) ./ (vecnorm(vec_shoulder_elbow,2,2).*vecnorm(vec_wrist_elbow,2,2)) );
        %throw out any discontinuities in angle changes
        ind = find(abs(diff(elbow_ang)) >= 360-30);
        for i = 1:length(ind)
            sign2pi = sign(elbow_ang(ind(i)+1)-elbow_ang(ind(i)));
            elbow_ang(ind(i):end) = elbow_ang(ind(i):end) - sign2pi*(360);
        end
        if min(elbow_ang) > 180 && max(elbow_ang) <= 360 %recenter the angle around 0
            elbow_ang = elbow_ang-360;
        elseif min(elbow_ang) < -360 && max(elbow_ang) < -180
            elbow_ang = elbow_ang+360;
        end
        
        
        %wrist rotation is the angle between the z axis of the wrist
        %tracker (in the standard coordinte frame) and the normal of the plane
        %for d = 1:size(BlockData{blk}.Grp(grp).mvmt(b,c).limbang.wrist,1)
        for d = 1:size(Yrsmp(:,:,6),1)
            %RotWrist = REuler(BlockData{blk}.Grp(grp).mvmt(b,c).limbang.wrist(d,1),BlockData{blk}.Grp(grp).mvmt(b,c).limbang.wrist(d,2),BlockData{blk}.Grp(grp).mvmt(b,c).limbang.wrist(d,3));
            %z_wrist(d,:) = RotWrist*[0 0 -1]';
            z_wrist(d,:) = Yrotangrsmp(:,:,d,5)*[0 0 1]';  %wrist rotation matrix times z vector at time d
        end
        wrist_ang = acosd( dot(n_armplane,z_wrist,2) ./ (vecnorm(n_armplane,2,2).*vecnorm(z_wrist,2,2)) );
        %throw out any discontinuities in angle changes
        ind = find(abs(diff(wrist_ang)) >= 360-30);
        for i = 1:length(ind)
            sign2pi = sign(wrist_ang(ind(i)+1)-wrist_ang(ind(i)));
            wrist_ang(ind(i):end) = wrist_ang(ind(i):end) - sign2pi*(360);
        end
        if min(wrist_ang) > 180 && max(wrist_ang) <= 360 %recenter the angle around 0
            wrist_ang = wrist_ang-360;
        elseif min(wrist_ang) < -360 && max(wrist_ang) < -180
            wrist_ang = wrist_ang+360;
        end
        
        Yarm = [ang_armplane_phi,ang_armplane_theta,ang_armplane_psi,elbow_ang,wrist_ang];
        
        
        
        %*********
        
        %Do a Procrustes comparison of the arm plane
        PDist.ArmPlane(a,c) = procrustesMultiScale(Xarm,Yarm,'reflect',0); %each row is a sample
        
%         figure(102)
%         subplot(6,7,c)
%         Xpca = Xarm - mean(Xarm);
%         [PCA_coef,PCA_score] = pca(Xpca,'centered',false);
%         Ypca = Yarm - mean(Yarm);
%         score = Ypca*PCA_coef;
%         
%         plot3(PCA_score(:,1),PCA_score(:,2),PCA_score(:,3),'k-');
%         hold on;
%         plot3(score(:,1),score(:,2),score(:,3),'b-');
%         hold off;
        
        
        %Do a Procrustes comparison of the joint angles
        PDist.JA(a,c) = procrustesMultiScale(Xjarsmp,Yjarsmp,'reflect',0); %each row is a sample
        
%         figure(103)
%         subplot(6,7,c)
%         Xpca = Xjarsmp - mean(Xjarsmp);
%         [PCA_coef,PCA_score] = pca(Xpca,'centered',false);
%         Ypca = Yjarsmp - mean(Yjarsmp);
%         score = Ypca*PCA_coef;
%         
%         plot3(PCA_score(:,1),PCA_score(:,2),PCA_score(:,3),'k-');
%         hold on;
%         plot3(score(:,1),score(:,2),score(:,3),'b-');
%         hold off;
        
        %Do a Procrustes comparison of the velocity matrix
        PDist.Vel(a,c) = procrustesMultiScale(Xvrsmpreshape,Yvrsmpreshape,'reflect',0); %each row is a sample
        
%         figure(104)
%         subplot(6,7,c)
%         Xpca = Xvrsmpreshape - mean(Xvrsmpreshape);
%         [PCA_coef,PCA_score] = pca(Xpca,'centered',false);
%         Ypca = Yvrsmpreshape - mean(Yvrsmpreshape);
%         score = Ypca*PCA_coef;
%         
%         plot3(PCA_score(:,1),PCA_score(:,2),PCA_score(:,3),'k-');
%         hold on;
%         plot3(score(:,1),score(:,2),score(:,3),'b-');
%         hold off;
        
        %Do a Procrustes comparison of the 3D velocity matrix
        PDist.Vel3D(a,c) = procrustesMultiScale(Xv3drsmp,Yv3drsmp,'reflect',0); %each row is a sample
        
%         figure(105)
%         subplot(6,7,c)
%         Xpca = Xv3drsmp - mean(Xv3drsmp);
%         [PCA_coef,PCA_score] = pca(Xpca,'centered',false);
%         Ypca = Yv3drsmp - mean(Yv3drsmp);
%         score = Ypca*PCA_coef;
%         
%         plot3(PCA_score(:,1),PCA_score(:,2),PCA_score(:,3),'k-');
%         hold on;
%         plot3(score(:,1),score(:,2),score(:,3),'b-');
%         hold off;
        
        %free up memory
        clear X* Y* ang* *ang vec* *hat n_armplane sew_axes sign2pi temp z_wrist;
        
    end

%     pause;
    
    
    %free up memory
    clear tmp;

    
end

%% plot 

fn = fieldnames(PDist);

for a = 1:length(fn)
    
    figure(a)
    clf;
    %histogram(PDist.(fn{a}));
    dat = PDist.(fn{a});
    histogram(nanmean(dat,2),[0:0.05:2]-0.025);
    
    title(fn{a});
    ylabel('Count');
    xlabel('Procrustes distance');
end



%% save


save(fullfile(fpath,'analyzedPdistData.mat'),'PDist','fnames');
