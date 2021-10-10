function outimage = SVDcompress(filename, Nretain)
    pic = imread(filename);
    % The image is loaded in as an integer format. Convert to floating point
    % so SVD can be taken.
    pic = double(pic);
    % The loaded data is in the range from 0-255. Convert to the range 0-1
    pic = pic / 255;
    
    outimage = zeros(height(pic), width(pic), 3);
    for i = 1:3
        % The loaded data is a 3-dimensional "matrix". The red, green and blue
        % components of the image are accessed with the last dimension.
        layer = pic(:,:,i);
        % Perform the SVD
        [U,S,V] = svd(layer);
        % The SVD is truncated by taking Nretain singular values 
        Uret = U(:,1:Nretain);
        Vret = V(:,1:Nretain);
        Sret = S(1:Nretain,1:Nretain);
        % Multiply the decomposition matrices back together
        reconstructed = Uret*Sret*Vret';
        % Save the reconstructed layer into the out image
        outimage(:,:,i) = reconstructed;
    end
end