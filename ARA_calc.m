







function ARA_calc(AllVolume,nblock)

tic;
[tp,nvox]=size(AllVolume);
Ara=zeros(1,nvox);
block_size=floor(nvox/nblock);
for bb=1:nblock
    left=1+(bb-1)*block_size;
    if bb==nblock
        right=nvox;
    else
        right=bb*block_size;
    end
    tmp_nvox=right-left+1;
    % local copy
    tmp_volume=AllVolume(:,left:right);
    % temporal smooth
    tmp_volume(2:tp-1,:)=0.25*tmp_volume(1:tp-2,:)+0.5*tmp_volume(2:tp-1,:)+0.25*tmp_volume(3:tp,:);
    tmp_volume(1,:)     =0.5*tmp_volume(1,:)+0.25*tmp_volume(2,:);
    tmp_volume(tp,:)    =0.25*tmp_volume(tp-1,:)+0.5*tmp_volume(tp,:);
    % current-former
    curr_fo=tmp_volume(2:tp-1,:)-tmp_volume(1:tp-2,:);
    curr_fo(curr_fo>0)= 1;
    curr_fo(curr_fo<0)=-1;
    % latter-current
    la_curr=tmp_volume(3:tp,:)-tmp_volume(2:tp-1,:);
    la_curr(la_curr>0)= 1;
    la_curr(la_curr<0)=-1;
    % peak-pit
    peak_pit=la_curr-curr_fo;
    clear la_curr curr_fo;
    % first point & last point
    fp=tmp_volume(2,:)-tmp_volume(1,:);
    fp_peak_idx=find(fp>0);
    fp_pit_idx =find(fp<0);
    fp(fp_peak_idx)=2; fp(fp_pit_idx)=-2;
    lp=tmp_volume(tp,:)-tmp_volume(tp-1,:);
    lp_peak_idx=find(lp>0);
    lp_pit_idx =find(lp<0);
    lp(lp_peak_idx)=-2; lp(lp_pit_idx)=2;
    peak_pit=[fp;peak_pit;lp];
    peak_idx=find(peak_pit==-2);
    pit_idx =find(peak_pit== 2);
    clear peak_pit;
    
    [peak_r,peak_c]=ind2sub([tp,tmp_nvox],peak_idx);
    npeak=hist(peak_c,1:tmp_nvox);
    [pit_r, pit_c] =ind2sub([tp,tmp_nvox],pit_idx);
    npit =hist(pit_c, 1:tmp_nvox);
    peak_data=tmp_volume(peak_idx);
    pit_data =tmp_volume(pit_idx);
    peak_theta=theta_ARA(peak_data,npeak);
    pit_theta =theta_ARA(pit_data, npit);
    Ara(left:right)=peak_theta./pit_theta;
end
Ara(isnan(Ara))=0; Ara(Ara==0)=1; Ara=log(Ara);
toc;

return
end