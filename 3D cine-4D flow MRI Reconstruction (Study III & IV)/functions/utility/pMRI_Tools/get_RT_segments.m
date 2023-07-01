function [max_inds] =get_RT_segments(flow,thresh)

% Segment the Heartbeats using flow information
dx = flow(1:end-1)-flow(2:end);
max_inds = find(dx>thresh*max(dx(:)));
for ii = 1:5
    inds_max_new = 1;
    ind2 = 1;
    while ind2<length(max_inds)+1
       if ind2 < length(max_inds)
           if abs(max_inds(ind2)-max_inds(ind2+1)) <3
               if dx(max_inds(ind2))>dx(max_inds(ind2+1))
                   max_inds_new(inds_max_new) = max_inds(ind2);
               else
                   max_inds_new(inds_max_new) = max_inds(ind2+1);
               end
               ind2 = ind2+2;
           else
               max_inds_new(inds_max_new) = max_inds(ind2);
               ind2 = ind2+1;
           end
           inds_max_new = inds_max_new+1;
       elseif ind2 == length(max_inds)
          if abs(max_inds(ind2)-max_inds(ind2-1)) <3
               if dx(max_inds(ind2))>dx(max_inds(ind2-1))
                   max_inds_new(inds_max_new) = max_inds(ind2);
               else
                   max_inds_new(inds_max_new) = max_inds(ind2-1);
               end
               ind2 = ind2+2;
           else
               max_inds_new(inds_max_new) = max_inds(ind2);
               ind2 = ind2+1;
          end
       end
    end
    max_inds = max_inds_new;
    clear max_inds_new
end
% crop_ind = round(150./(1000*2*TR*lpf));
% max_inds = max(max_inds-crop_ind(ind),1);