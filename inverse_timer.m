
tic;
counter = 0;
while toc < 1
    a = rand(100, 100);
    b = inv(a);
    counter = counter + 1;
end

disp(counter);