function im = imImport(filename)
    im = imread(filename);
    im = uint16(im(:,:,2))*8 + uint16(im(:,:,3))/8;
    im = im;
end