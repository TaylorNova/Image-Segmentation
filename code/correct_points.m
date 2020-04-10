function c_point = correct_points(point,height,width)
    if point(1) < 1
        xlabel = 1;
    elseif point(1) > height
        xlabel = height;
    else
        xlabel = point(1);
    end
    if point(2) < 1
        ylabel = 1;
    elseif point(2) > width
        ylabel = width;
    else
        ylabel = point(2);
    end
    c_point = [xlabel,ylabel];
end