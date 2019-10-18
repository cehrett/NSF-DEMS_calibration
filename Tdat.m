% Rescale inputs and standardize outputs from simulator

function tdat = Tdat(dat,num_out)
% dat is a matrix containing the input columns followed by the
% response columns. num_out is the number of trailing columns that
% are response columns. dat contains exactly one control parameter, as its
% first column.

% Get number of total columns and rows
cols = size(dat,2);
rows = size(dat,1);

% Get input columns
dat_in = dat(:,1:(cols-num_out));

% Get output columns
dat_out = dat(:,(cols-num_out+1):end);

% Add dummy var input cols
%dum_var = repelem(1:num_out,rows)';
dum_mat = [eye(num_out-1) ; zeros(1,num_out-1)];
if num_out > 1
    dum_mat = repelem(dum_mat,rows,1);
    dat_in = [ dum_mat repmat(dat_in,num_out,1) ];
end

% Rescale each input column
tdat_in = zeros(size(dat_in)); % this will store rescaled input
in_cols = size(dat_in,2);
for ii = 1:in_cols
    input_mins(ii)   = min(dat_in(:,ii));
    input_ranges(ii) = range(dat_in(:,ii));
    tdat_in(:,ii) = (dat_in(:,ii) - input_mins(ii))/input_ranges(ii);
end

% % Standardize each output column for each temp
% % This version is good for outputs collected on a grid over the control
% inputs, at a discrete number of control input locations at each of which
% many observations were made
% tdat_out = zeros(size(dat_out)); % this will store standardized output
% for ii = 1:num_out
%     settings = unique(dat(:,1)); % We want the mean for each setting
%                                  % of the first parameter (temperature)
%     for jj = 1:length(settings)
%         dat_set = dat_out(dat(:,1)==settings(jj),:);
%         % dat_set is the output at the current setting (temp)
%         output_means(ii,jj) = mean(dat_set(:,ii));
%         output_sds(ii,jj)   = std(dat_set(:,ii));
%         tdat_out(dat(:,1)==settings(jj),ii) = (dat_set(:,ii) - ...
%             output_means(ii,jj))/output_sds(ii,jj);
%     end
% end

% Alt: standardize each output column across all temps
tdat_out = zeros(size(dat_out)); % this will store standardized output
for ii = 1:num_out
    output_means(ii,:) = mean(dat_out(:,ii)) ;
    output_sds(ii,:) = std(dat_out(:,ii)) ;
    tdat_out(:,ii) = (dat_out(:,ii) - output_means(ii,:)) / ...
        output_sds(ii,:) ;
end

% Get single output column
tdat_out = tdat_out(:);

% Pack up and leave
tdat.input          = tdat_in;
tdat.output         = tdat_out;
tdat.input_mins     = input_mins;
tdat.input_ranges   = input_ranges;
tdat.output_means   = output_means;
tdat.output_sds     = output_sds;
tdat.cntrl_settings = settings;

end
