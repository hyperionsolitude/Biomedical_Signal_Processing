close all
clear
clc
%% Question 1
%% a)

img1=im2double(imread("Fig3.08(a).jpg"));
power_law=@(c,r,gamma) real(c.*(r.^gamma));
out_img1_06=power_law(1,img1,0.6);
out_img1_04=power_law(1,img1,0.4);
out_img1_03=power_law(1,img1,0.3);

figure;
subplot(221);
imshow(img1);
title("Original Image");
subplot(222);
imshow(out_img1_06);
title("Power Law c=1,gamma=0.6");
subplot(223);
imshow(out_img1_04);
title("Power Law c=1,gamma=0.4");
subplot(224);
imshow(out_img1_03);
title("Power Law c=1,gamma=0.3");

%% b)

%% i

a=im2double(imread("Fig3.46(a).jpg"));
figure;
imshow(a);
title("Original Image (a)");


%% ii

H_laplacian=fspecial('laplacian');
b=imfilter(a,-H_laplacian);
figure;
imshow(b);
title("Laplacian Filtering of a given Image");

%% iii

c=a+b;
figure;
imshow(c);
title("Sharpened Image");

%% iv

H_sobel=fspecial('sobel');
d=imfilter(a,H_sobel);
figure;
imshow(d);
title("Sobel Filtering of a given Image");

%% v

H_averaging=fspecial('average',[5,5]);
e=imfilter(d,H_averaging);
figure;
imshow(e);
title("Averaging Filtering of a given Image");

%% vi

f=c.*e;
figure;
imshow(f);
title("Mask of a given Image by product of c and e");

%% vii

g=a+f;
figure;
imshow(g);
title("Sharpened Image by summation of a and f");

%% viii

h=power_law(1,g,0.7);
figure;
subplot(251)
imshow(power_law(1,g,0.1));
title("Final Image by Power Law c=1,gamma=0.1");
subplot(252)
imshow(power_law(1,g,0.2));
title("Final Image by Power Law c=1,gamma=0.2");
subplot(253)
imshow(power_law(1,g,0.3));
title("Final Image by Power Law c=1,gamma=0.3");
subplot(254)
imshow(power_law(1,g,0.4));
title("Final Image by Power Law c=1,gamma=0.4");
subplot(255)
imshow(power_law(1,g,0.5));
title("Final Image by Power Law c=1,gamma=0.5");
subplot(256)
imshow(power_law(1,g,0.6));
title("Final Image by Power Law c=1,gamma=0.6");
subplot(257)
imshow(power_law(1,g,0.7));
title("Final Image by Power Law c=1,gamma=0.7");
subplot(258)
imshow(power_law(1,g,0.8));
title("Final Image by Power Law c=1,gamma=0.8");
subplot(259)
imshow(power_law(1,g,0.9));
title("Final Image by Power Law c=1,gamma=0.9");
subplot(2,5,10)
imshow(power_law(1,g,1));
title("Final Image by Power Law c=1,gamma=1");

%% ix

figure;
subplot(131);
imshow(g);
title("Sharpened Image g");
subplot(132);
imshow(h);
title("Power Law Resulting Image h");
subplot(133);
imshow(a);
title("Original Image a");

%% Question 2
I=im2double(imread("Fig0462(a)(PET_image).tif"));
% Get ln of Image

I = log(1 + I);

% Get M, N for FFT

[M,N] = size(I);

% High Pass Filter

sigma = 20;
const=1;
Gamma_h=3.0;
Gamma_l=0.4;

p=ceil(M/2);
q=ceil(N/2);
H=zeros(M,N);
for i=1:M
for j=1:N
distance=sqrt((i-p)^2+(j-q)^2);
H(i,j)=(Gamma_h-Gamma_l)*(1-exp((-const*(distance)^2)/((sigma^2))))+Gamma_l;
end
end

% Centering Filter Response

H = fftshift(H);

% FFT

If = fft2(I, M, N);

% Filtering Then Taking Inverse FFT

Iout = real(ifft2(H.*If));
Iout = Iout(1:size(I,1),1:size(I,2));

% Taking Inverse Logarithm

Ihmf = exp(Iout) - 1;

% display the images
figure;
imshowpair(I, Ihmf, 'montage')
title("Original image vs Enchanced image");