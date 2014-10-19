


    tp=size(seed,1);
    % temporal smooth
%     seed(2:tp-1,:)=0.25*seed(1:tp-2,:)+0.5*seed(2:tp-1,:)+0.25*seed(3:tp,:);
%     seed(1,:)     =0.5*seed(1,:)+0.25*seed(2,:);
%     seed(tp,:)    =0.25*seed(tp-1,:)+0.5*seed(tp,:);
    %- 2:end-1 point
    % current-former
    curr_fo=seed(2:tp-1,:)-seed(1:tp-2,:);
    curr_fo(curr_fo>0)= 1;
    curr_fo(curr_fo<0)=-1;
    % latter-current
    la_curr=seed(3:tp,:)-seed(2:tp-1,:);
    la_curr(la_curr>0)= 1;
    la_curr(la_curr<0)=-1;
    % peak-pit
    peak_pit=la_curr-curr_fo;
    %- first point & last point
    fp=seed(2,:)-seed(1,:);
    fp_peak_idx=find(fp>0);
    fp_pit_idx =find(fp<0);
    fp(fp_peak_idx)=2; fp(fp_pit_idx)=-2;
    lp=seed(tp,:)-seed(tp-1,:);
    lp_peak_idx=find(lp>0);
    lp_pit_idx =find(lp<0);
    lp(lp_peak_idx)=-2; lp(lp_pit_idx)=2;
    peak_pit=[fp;peak_pit;lp];
    peak_idx=find(peak_pit==-2);
    pit_idx =find(peak_pit== 2);
 
    
    
    plot(1:tp,seed,'k'); hold on; peak_plot=zeros(tp,1);plot(peak_idx,seed(peak_idx),'ro'); plot(pit_idx,seed(pit_idx),'bo'); hold off;
    
    [peak_r,peak_c]=ind2sub([tp,1],peak_idx);
    npeak=hist(peak_c,1:1);
    [pit_r, pit_c] =ind2sub([tp,1],pit_idx);
    npit =hist(pit_c, 1:1);
    peak_data=seed(peak_idx);
    pit_data =seed(pit_idx);
    peak_theta=theta_ARA(peak_data,npeak);
    pit_theta =theta_ARA(pit_data, npit);
    seed_Ara=peak_theta./pit_theta;
    seed_Ara(isnan(seed_Ara))=0;