

data_files=spm_select(Inf,'any','Select Data Files');
mask_files=spm_select(Inf,'any','Sleect Mask Files');

data_hdr=spm_vol(data_files(1,:));
% mask intersect
mask_dat=zeros(data_hdr.dim);
for ii=1:length(mask_files)
    tmp_mask_fname=mask_files(ii,:);
    tmp_mask_hdr=spm_vol(tmp_mask_file);
    if ~all(data_hdr.mat==tmp_mask_hdr.mat)
        fprintf('mask %d dimension not match, resizing now\n',ii);
        tmp_arg=struct('mean',false,'interp' ,0,'which',1, 'prefix','tmp_' );
        spm_reslice([data_hdr tmp_mask_hdr],tmp_arg); %会在N头文件所在目录下产生一个文件
        [tmp_mask_path,tmp_mask_fname]=fileparts(tmp_mask_fname);
        tmp_mask_fname=[tmp_mask_path,filesep,tmp_arg.prefix,tmp_mask_fname];
        tmp_mask_dat=logical(spm_read_vols(spm_vol(tmp_mask_fname)));
        delete(tmp_mask_fname);
        clear tmp_mask_fname;
    else
        tmp_mask_dat=logical(spm_read_vols(tmp_mask_hdr));
    end
    mask_dat=mask_dat+tmp_mask_dat;
end
mask_idx=find(logical(mask_dat)); clear mask_dat; 

%# load data and compute kernel Phi
nvox=sum(mask_dat(:));
nfile=length(data_files);
% memory limits
def=prt_get_defaults('fs');
mem=def.mem_limit;
block_size=floor(mem/(8*3)/nfile);
n_block= ceil(nvox/block_size);
if block_size<nvox
    beep; fprintf('out of memory\n'); pause;
end
h = waitbar(0,'Please wait while preparing feature set');
feature_mat=zeros(nfile,nvox);
Phi=zeros(nfile);
bstart=1; bend=min(nvox,block_size);
for ii=1:n_block
    disp ([' > preparing block: ', num2str(b),' of ',num2str(n_block),' ...']);
    tmp_vox_range=bstart:bend;
    tmp_vox_num=length(tmp_vox_range);
    tmp_vox_idx=mask_idx(tmp_vox_range);
    kern_vols=zeros(nfile,tmp_vox_num);
    for jj=1:nfile
        tmp_data_hdr=spm_vol(data_files(jj,:));
        if jj>1&&~all(tmp_data_hdr.mat==data_hdr.mat)
           error('data dimension not match');
        end
        tmp_dat=spm_read_vols(tmp_data_hdr);
        tmp_dat(isnan(tmp_dat))=0; tmp_dat=tmp_dat(tmp_vox_idx);
        kern_vols(ii,:)=tmp_dat;
    end
    % size limit setup such that no more than ~1Gb of mem is required:
    % 1Gb/3(nr of matrices)/8(double)= ~40e6 -> sqrt -> 6e3 element
    if size(kern_vols,1)>6e3
        for ic=1:size(kern_vols,1)
            Phi(ic,:)=Phi(ic,:)+kern_vols*kern_vols(ic,:)';
        end
    else
        Phi=Phi+kern_vols*kern_vols';
    end
    bstart = bend+1; bend = min(bstart+block_size-1,n_vox);
    clear kern_vols;
end
close(h);    
  
