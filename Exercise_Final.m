%% Final assignment
clear variables
close all

video = VideoReader('MAH01462.MP4');
videoWidth = video.Width;
videoHeight = video.Height;

while hasFrame(video)
    frame = readFrame(video, 'native');
    imshow(frame)
end

%getpts();