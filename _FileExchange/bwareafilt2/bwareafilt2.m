% Function BWAREAFILT2, replacement for older versions of matlab
% Created by Robin Bruggink / robinbruggink@gmail.com
% Input BW (2D or 3D image in BW)
%       range (The range of area size to search within), [] will allow
%       maximum range
%       p (Number to keep)
%       b (Largest of Smallest)
function result = bwareafilt2(BW, range, p, d)

    % Check input
    if(nargin ~= 4)
        error('Use 4 input arguments, bwareafilt2(BW (bw image or volume),range (array of with the lower and larger range),p (number to keep),d (smallest of largest to keep)');
    end
   
    % Convert the input to logical
    BW = BW > 0.5;
    
    % Get the regionprops
    props = regionprops(BW,'area','pixelidxlist');
    
    % Remove the areas outside the range
    if(~isempty(range))
        idx = [find([props.Area] < range(1)) find([props.Area] > range(2))];
        props(idx) = [];
    end
    
    % If the regionprops count is less than the input, give a warning
    if(length(props) < p)
        p = length(props);
        warning('There were less areas found than provided in the input');
    end
    
    % Sort the regionprops to the right order (smallest or largest)
    switch d
        case 'smallest'
            [~, ind]=sort([props.Area],2,'ascend');
        case 'largest'
            [~, ind]=sort([props.Area],2,'descend');
        otherwise
            [~, ind]=sort([props.Area],2,'ascend');
    end
    % Sort
    props=props(ind);
    
    assignin('base','prop',props);
    
    % Retrieve the n largest of smallest objects.
    result = zeros(size(BW));
    for i = 1:p
       result(props(i).PixelIdxList) = 1; 
    end
    
    % Convert the result to logical
    result = result > 0.5;
end