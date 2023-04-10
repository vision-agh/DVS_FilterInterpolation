classdef EventFilterInterpolation < handle
% EventFilterSimpleBlocksAdvancedInterpolation(frame_size, filter_length, scale, update_factor, filtered_ts, interpolation_method)
% Interpolation methods:
% 0 - bilinear
% 1 - bilinear with interval weights
% 2 - max
% 3 - distance

    properties
        FrameSize
        FilterLength
        Scale
        UpdateFactor
        TimestampMap
        IntervalMap
        ActiveMap
        EventsTrue
        EventsFalse
        ValidEvents
        InvalidEvents
        CurrentTs
        ImageData
        InterpolationMethod
    end
    
    methods
        function obj = EventFilterSimpleBlocksAdvancedInterpolation(frame_size, filter_length, scale, update_factor, filtered_ts, interpolation_method)
            obj.FrameSize = frame_size;
            obj.FilterLength = filter_length;
            obj.Scale = scale;
            obj.UpdateFactor = update_factor;
            obj.TimestampMap = zeros(floor(frame_size(1)/scale), floor(frame_size(2)/scale));
            obj.IntervalMap = 1e4 * ones(floor(frame_size(1)/scale), floor(frame_size(2)/scale));
            obj.ActiveMap = zeros(floor(frame_size(1)/scale), floor(frame_size(2)/scale));
            obj.ValidEvents = 0;
            obj.InvalidEvents = 0;
            obj.CurrentTs = 0;
            obj.FilteredTs = filtered_ts;
            obj.InterpolationMethod = interpolation_method;
        end

        function events = getEventsTrue(obj)
            events = obj.EventsTrue;
        end

        function events = getEventsFalse(obj)
            events = obj.EventsFalse;
        end
        
        function processEvents(obj, events)
            eventsBin = false(1, size(events, 1));
            for iter = 1:size(events, 1)
                event = events(iter, :);
                eventsBin(iter) = obj.filterEvent(event);
                obj.updateFeatures(event);
                obj.CurrentTs = event(4);
            end
            obj.updateEmpty()
            obj.EventsTrue = events(eventsBin, :);
            obj.EventsFalse = events(~eventsBin, :);
            obj.ValidEvents = obj.ValidEvents + size(obj.EventsTrue, 1);
            obj.InvalidEvents = obj.InvalidEvents + size(obj.EventsFalse, 1);
        end

        function correct = filterEvent(obj, event)
            switch obj.InterpolationMethod
                case 0
                    thr_ts = obj.bilinear_interpolation(event(2), event(1));
                case 1
                    thr_ts = obj.bilinear_interval_weights_interpolation(event(2), event(1));
                case 2
                    thr_ts = obj.max_interpolation(event(2), event(1));
                case 3
                    thr_ts = obj.distance_interpolation_init(event(2), event(1));
            end
            diff_ts = event(4) - thr_ts;
            correct = diff_ts < obj.FilterLength;
        end

        function updateFeatures(obj, event)
            obj.updateInterval(event);
            obj.updateFilteredTimestamp(event);
            obj.updateActive(event);
        end

        function updateInterval(obj, event)
            cellY = floor(event(2)/obj.Scale)+1;
            cellX = floor(event(1)/obj.Scale)+1;
            cell_time = obj.TimestampMap(cellY, cellX);
            cell_interval = obj.IntervalMap(cellY, cellX);
            time_diff = event(4) - cell_time;
            new_interval = cell_interval * (1 - obj.UpdateFactor) + time_diff * obj.UpdateFactor;
            obj.IntervalMap(cellY, cellX) = new_interval;
        end

        function updateFilteredTimestamp(obj, event)
            cellY = floor(event(2)/obj.Scale)+1;
            cellX = floor(event(1)/obj.Scale)+1;
            cell_filtered_time = obj.TimestampMap(cellY, cellX);
            new_filtered_time = cell_filtered_time * (1 - obj.UpdateFactor) + event(4) * obj.UpdateFactor;
            obj.TimestampMap(cellY, cellX) = new_filtered_time;
        end

        function updateActive(obj, event)
            cellY = floor(event(2)/obj.Scale)+1;
            cellX = floor(event(1)/obj.Scale)+1;
            obj.ActiveMap(cellY, cellX) = 1;
        end

        function updateEmpty(obj)
            scaled_size = size(obj.ActiveMap);
            for y = 1:scaled_size(1)
                for x = 1:scaled_size(2)
                    if ~obj.ActiveMap(y, x)
                        % Update interval
                        cell_time = obj.TimestampMap(y, x);
                        cell_interval = obj.IntervalMap(y, x);
                        time_diff = obj.CurrentTs - cell_time;
                        new_interval = cell_interval * (1 - obj.UpdateFactor) + time_diff * obj.UpdateFactor;
                        obj.IntervalMap(y, x) = new_interval;
                        % Update filtered timestamp
                        new_filtered_time = cell_time * (1 - obj.UpdateFactor) + obj.CurrentTs * obj.UpdateFactor;
                        obj.TimestampMap(y, x) = new_filtered_time;
                    end
                end
            end
            obj.ActiveMap(:, :, 6) = 0;
        end
        
        function displayFilteredTimestamp(obj)
            diff_map = obj.CurrentTs - obj.TimestampMap;
            exp_map = 255 * exp(-diff_map / 100);
            exp_map = imresize(exp_map, obj.FrameSize, "nearest");
            figure(2);
            imshow(exp_map);
            title('Filtered timestamp');
        end
        
        function exp_map = returnFilteredTimestamp(obj)
            diff_map = obj.CurrentTs - obj.TimestampMap;
            exp_map = 255 * exp(-diff_map / 100);
            exp_map = imresize(exp_map, obj.FrameSize, "nearest");
        end
        
        function displayInterval(obj)
            exp_map = 255 * exp(-obj.IntervalMap / 20);
            exp_map = imresize(exp_map, obj.FrameSize, "nearest");
            figure(3)
            imshow(exp_map);
            title('Event intervals');
        end
        
        function exp_map = returnInterval(obj)
            exp_map = 255 * exp(-obj.IntervalMap / 20);
            exp_map = imresize(exp_map, obj.FrameSize, "nearest");
        end

        function displayFeatures(obj)
            image_rgb = zeros(obj.FrameSize(1), obj.FrameSize(2), 3);
            image_rgb(:, :, 1) = obj.returnFilteredTimestamp();
            image_rgb(:, :, 2) = obj.returnInterval();
            figure(4);
            imshow(image_rgb);
            title('Features');
        end

        function displayFilter(obj)
            image = uint8(ones(obj.FrameSize(1), obj.FrameSize(2), 3)) * 255;
            if ~isempty(obj.EventsFalse)
                for ind = 1:size(obj.EventsFalse, 1)
                    image(obj.EventsFalse(ind, 2)+1, obj.EventsFalse(ind, 1)+1, :) = [255, 0, 0];
                end
            end
            if ~isempty(obj.EventsTrue)
                for ind = 1:size(obj.EventsTrue, 1)
                    image(obj.EventsTrue(ind, 2)+1, obj.EventsTrue(ind, 1)+1, :) = [0, 255, 0];
                end
            end
            figure(5);
            imshow(image);
            obj.ImageData = image;
        end

        function saveImage(obj, fileName)
            image = uint8(ones(obj.FrameSize(1), obj.FrameSize(2), 3)) * 255;
            if ~isempty(obj.EventsFalse)
                for ind = 1:size(obj.EventsFalse, 1)
                    image(obj.EventsFalse(ind, 2)+1, obj.EventsFalse(ind, 1)+1, :) = [255, 0, 0];
                end
            end
            if ~isempty(obj.EventsTrue)
                for ind = 1:size(obj.EventsTrue, 1)
                    image(obj.EventsTrue(ind, 2)+1, obj.EventsTrue(ind, 1)+1, :) = [0, 255, 0];
                end
            end
            imwrite(image, fileName);
        end

        function displaySaveFilter(obj, fileName)
            image = uint8(ones(obj.FrameSize(1), obj.FrameSize(2), 3)) * 255;
            if ~isempty(obj.EventsFalse)
                for ind = 1:size(obj.EventsFalse, 1)
                    image(obj.EventsFalse(ind, 2)+1, obj.EventsFalse(ind, 1)+1, :) = [255, 0, 0];
                end
            end
            if ~isempty(obj.EventsTrue)
                for ind = 1:size(obj.EventsTrue, 1)
                    image(obj.EventsTrue(ind, 2)+1, obj.EventsTrue(ind, 1)+1, :) = [0, 255, 0];
                end
            end
            figure(5);
            imshow(image);
            obj.ImageData = image;
            imwrite(image, fileName);
        end
        
        function thr = bilinear_interval_weights_interpolation(obj, y, x)
            y_cell = floor(y / obj.Scale)+1;
            x_cell = floor(x / obj.Scale)+1;

            w = obj.Scale;

            if (y < obj.Scale / 2) || (y >= obj.FrameSize(1) - obj.Scale / 2)
                if (x < obj.Scale / 2) || (x >= obj.FrameSize(2) - obj.Scale / 2)

                    ts_filtered_cell = obj.TimestampMap(y_cell, x_cell);
                    thr = ts_filtered_cell;

                else
                    i_floor = floor((x - w / 2) / w);
                    dx1 = x - w / 2 - i_floor * w + 0.5;
                    dx2 = (i_floor + 1) * w - (x - w / 2) - 0.5;

                    interval_11 = obj.IntervalMap(y_cell, i_floor + 1);
                    interval_12 = obj.IntervalMap(y_cell, i_floor + 2);

                    coef_11 = dx2 * interval_12;
                    coef_12 = dx1 * interval_11;

                    sum_top = coef_11 + coef_12;

                    ts_filtered_cell_11 = obj.TimestampMap(y_cell, i_floor + 1);
                    ts_filtered_cell_12 = obj.TimestampMap(y_cell, i_floor + 2);

                    thr = ts_filtered_cell_11 * coef_11 / sum_top + ts_filtered_cell_12 * coef_12 / sum_top;
                    
                end
            else
                if (x < obj.Scale / 2) || (x >= obj.FrameSize(2) - obj.Scale / 2)
                    j_floor = floor((y - w / 2) / w);
                    dy1 = y - w / 2 - j_floor * w + 0.5;
                    dy2 = (j_floor + 1) * w - (y - w / 2) - 0.5;

                    interval_11 = obj.IntervalMap(j_floor + 1, x_cell);
                    interval_21 = obj.IntervalMap(j_floor + 2, x_cell);

                    coef_top = dy2 * interval_21;
                    coef_bot = dy1 * interval_11;

                    sum_all = coef_top + coef_bot;

                    ts_filtered_cell_11 = obj.TimestampMap(j_floor + 1, x_cell);
                    ts_filtered_cell_21 = obj.TimestampMap(j_floor + 2, x_cell);

                    thr = ts_filtered_cell_11 * coef_top / sum_all + ts_filtered_cell_21 * coef_bot / sum_all;

                else
                    i_floor = floor((x - w / 2) / w);
                    j_floor = floor((y - w / 2) / w);
                    dx1 = x - w / 2 - i_floor * w + 0.5;
                    dx2 = (i_floor + 1) * w - (x - w / 2) - 0.5;
                    dy1 = y - w / 2 - j_floor * w + 0.5;
                    dy2 = (j_floor + 1) * w - (y - w / 2) - 0.5;

                    interval_11 = obj.IntervalMap(j_floor + 1, i_floor + 1);
                    interval_12 = obj.IntervalMap(j_floor + 1, i_floor + 2);
                    interval_21 = obj.IntervalMap(j_floor + 2, i_floor + 1);
                    interval_22 = obj.IntervalMap(j_floor + 2, i_floor + 2);

                    interval_top = interval_11 * interval_12;
                    interval_bot = interval_21 * interval_22;

                    coef_11 = dx2 * interval_12;
                    coef_12 = dx1 * interval_11;
                    coef_21 = dx2 * interval_22;
                    coef_22 = dx1 * interval_21;

                    sum_top = coef_11 + coef_12;
                    sum_bot = coef_21 + coef_22;

                    coef_top = dy2 * interval_bot;
                    coef_bot = dy1 * interval_top;
                    sum_all = coef_top + coef_bot;
                    
                    ts_filtered_cell_11 = obj.IntervalMap(j_floor + 1, i_floor + 1);
                    ts_filtered_cell_12 = obj.IntervalMap(j_floor + 1, i_floor + 2);
                    ts_filtered_cell_21 = obj.IntervalMap(j_floor + 2, i_floor + 1);
                    ts_filtered_cell_22 = obj.IntervalMap(j_floor + 2, i_floor + 2);

                    thr1 = ts_filtered_cell_11 * coef_11 / sum_top + ts_filtered_cell_12 * coef_12 / sum_top;
                    thr2 = ts_filtered_cell_21 * coef_21 / sum_bot + ts_filtered_cell_22 * coef_22 / sum_bot;
                    thr = thr1 * coef_top / sum_all + thr2 * coef_bot / sum_all;

                end
            end
        end

        function thr = bilinear_interpolation(obj, y, x)
            y_cell = floor(y / obj.Scale)+1;
            x_cell = floor(x / obj.Scale)+1;

            w = obj.Scale;

            if (y < obj.Scale / 2) || (y >= obj.FrameSize(1) - obj.Scale / 2)
                if (x < obj.Scale / 2) || (x >= obj.FrameSize(2) - obj.Scale / 2)

                    ts_filtered_cell = obj.TimestampMap(y_cell, x_cell);
                    thr = ts_filtered_cell;
                else
                    i_floor = floor((x - w / 2) / w);
                    dx1 = x - w / 2 - i_floor * w + 0.5;
                    dx2 = (i_floor + 1) * w - (x - w / 2) - 0.5;

                    ts_filtered_cell_11 = obj.TimestampMap(y_cell, i_floor + 1);
                    ts_filtered_cell_12 = obj.TimestampMap(y_cell, i_floor + 2);

                    thr = ts_filtered_cell_11 * (dx2 / w) + ts_filtered_cell_12 * (dx1 / w);

                end
            else
                if (x < obj.Scale / 2) || (x >= obj.FrameSize(2) - obj.Scale / 2)
                    j_floor = floor((y - w / 2) / w);
                    dy1 = y - w / 2 - j_floor * w + 0.5;
                    dy2 = (j_floor + 1) * w - (y - w / 2) - 0.5;

                    ts_filtered_cell_11 = obj.TimestampMap(j_floor + 1, x_cell);
                    ts_filtered_cell_21 = obj.TimestampMap(j_floor + 2, x_cell);

                    thr = ts_filtered_cell_11 * (dy2 / w) + ts_filtered_cell_21 * (dy1 / w);

                else
                    i_floor = floor((x - w / 2) / w);
                    j_floor = floor((y - w / 2) / w);
                    dx1 = x - w / 2 - i_floor * w + 0.5;
                    dx2 = (i_floor + 1) * w - (x - w / 2) - 0.5;
                    dy1 = y - w / 2 - j_floor * w + 0.5;
                    dy2 = (j_floor + 1) * w - (y - w / 2) - 0.5;

                    ts_filtered_cell_11 = obj.TimestampMap(j_floor + 1, i_floor + 1);
                    ts_filtered_cell_12 = obj.TimestampMap(j_floor + 1, i_floor + 2);
                    ts_filtered_cell_21 = obj.TimestampMap(j_floor + 2, i_floor + 1);
                    ts_filtered_cell_22 = obj.TimestampMap(j_floor + 2, i_floor + 2);

                    thr1 = ts_filtered_cell_11 * (dx2 / w) + ts_filtered_cell_12 * (dx1 / w);
                    thr2 = ts_filtered_cell_21 * (dx2 / w) + ts_filtered_cell_22 * (dx1 / w);
                    thr = thr1 * (dy2 / w) + thr2 * (dy1 / w);

                end
            end
        end

        function thr = max_interpolation(obj, y, x)
            y_cell = floor(y / obj.Scale)+1;
            x_cell = floor(x / obj.Scale)+1;

            w = obj.Scale;

            if (y < obj.Scale / 2) || (y >= obj.FrameSize(1) - obj.Scale / 2)
                if (x < obj.Scale / 2) || (x >= obj.FrameSize(2) - obj.Scale / 2)

                    ts_filtered_cell = obj.TimestampMap(y_cell, x_cell);
                    thr = ts_filtered_cell;
                else
                    i_floor = floor((x - w / 2) / w);

                    ts_filtered_cell_11 = obj.TimestampMap(y_cell, i_floor + 1);
                    ts_filtered_cell_12 = obj.TimestampMap(y_cell, i_floor + 2);

                    thr = max(ts_filtered_cell_11, ts_filtered_cell_12);

                end
            else
                if (x < obj.Scale / 2) || (x >= obj.FrameSize(2) - obj.Scale / 2)
                    j_floor = floor((y - w / 2) / w);

                    ts_filtered_cell_11 = obj.TimestampMap(j_floor + 1, x_cell);
                    ts_filtered_cell_21 = obj.TimestampMap(j_floor + 2, x_cell);

                    thr = max(ts_filtered_cell_11, ts_filtered_cell_21);

                else
                    i_floor = floor((x - w / 2) / w);
                    j_floor = floor((y - w / 2) / w);
                    
                    ts_filtered_cell_11 = obj.TimestampMap(j_floor + 1, i_floor + 1);
                    ts_filtered_cell_12 = obj.TimestampMap(j_floor + 1, i_floor + 2);
                    ts_filtered_cell_21 = obj.TimestampMap(j_floor + 2, i_floor + 1);
                    ts_filtered_cell_22 = obj.TimestampMap(j_floor + 2, i_floor + 2);

                    thr = max([ts_filtered_cell_11, ts_filtered_cell_12, ts_filtered_cell_21, ts_filtered_cell_22]);

                end
            end
        end

        function thr = distance_interpolation(obj, y, x)
            y_cell = floor(y / obj.Scale)+1;
            x_cell = floor(x / obj.Scale)+1;

            w = obj.Scale;

            if (y < obj.Scale / 2) || (y >= obj.FrameSize(1) - obj.Scale / 2)
                if (x < obj.Scale / 2) || (x >= obj.FrameSize(2) - obj.Scale / 2)

                    ts_filtered_cell = obj.TimestampMap(y_cell, x_cell);
                    thr = ts_filtered_cell;
                else
                    i_floor = floor((x - w / 2) / w);
                    dx1 = x - w / 2 - i_floor * w + 0.5;
                    dx2 = (i_floor + 1) * w - (x - w / 2) - 0.5;

                    if (y < obj.Scale / 2)
                    dy =  w / 2 - y - 0.5;
                    else
                    j_floor = floor((y - w / 2) / w);
                    dy = y - w / 2 - j_floor * w + 0.5;
                    end

                    dist_1 = sqrt(dx1^2+dy^2);
                    dist_2 = sqrt(dx2^2+dy^2);

                    interval_1 = obj.IntervalMap(y_cell, i_floor + 1);
                    interval_2 = obj.IntervalMap(y_cell, i_floor + 2);

                    C_1 = dist_2 * interval_2;
                    C_2 = dist_1 * interval_1;
        
                    coef_sum = C_1 + C_2;

                    ts_filtered_cell_11 = obj.TimestampMap(y_cell, i_floor + 1);
                    ts_filtered_cell_12 = obj.TimestampMap(y_cell, i_floor + 2);

                    thr = (ts_filtered_cell_11 * C_1 + ts_filtered_cell_12 * C_2) / coef_sum;

                end
            else
                if (x < obj.Scale / 2) || (x >= obj.FrameSize(2) - obj.Scale / 2)
                    if (x < obj.Scale / 2)
                        dx = w / 2 - x - 0.5;
                    else
                        i_floor = floor((x - w / 2) / w);
                        dx = x - w / 2 - i_floor * w + 0.5;
                    end

                    j_floor = floor((y - w / 2) / w);
                    dy1 = y - w / 2 - j_floor * w + 0.5;
                    dy2 = (j_floor + 1) * w - (y - w / 2) - 0.5;

                    dist_1 = sqrt(dx^2+dy1^2);
                    dist_2 = sqrt(dx^2+dy2^2);

                    interval_1 = obj.IntervalMap(j_floor + 1, x_cell);
                    interval_2 = obj.IntervalMap(j_floor + 2, x_cell);

                    C_1 = dist_2 * interval_2;
                    C_2 = dist_1 * interval_1;

                    coef_sum = C_1 + C_2;

                    ts_filtered_cell_11 = obj.TimestampMap(j_floor + 1, x_cell);
                    ts_filtered_cell_21 = obj.TimestampMap(j_floor + 2, x_cell);

                    thr = (ts_filtered_cell_11 * C_1 + ts_filtered_cell_21 * C_2) / coef_sum;

                else
                    i_floor = floor((x - w / 2) / w);
                    j_floor = floor((y - w / 2) / w);
                    dx1 = x - w / 2 - i_floor * w + 0.5;
                    dx2 = (i_floor + 1) * w - (x - w / 2) - 0.5;
                    dy1 = y - w / 2 - j_floor * w + 0.5;
                    dy2 = (j_floor + 1) * w - (y - w / 2) - 0.5;

                    dist_11 = sqrt(dx1^2+dy1^2);
                    dist_12 = sqrt(dx2^2+dy1^2);
                    dist_21 = sqrt(dx1^2+dy2^2);
                    dist_22 = sqrt(dx2^2+dy2^2);

                    interval_11 = obj.IntervalMap(j_floor + 1, i_floor + 1);
                    interval_12 = obj.IntervalMap(j_floor + 1, i_floor + 2);
                    interval_21 = obj.IntervalMap(j_floor + 2, i_floor + 1);
                    interval_22 = obj.IntervalMap(j_floor + 2, i_floor + 2);

                    coef_11 = dist_11 * interval_11;
                    coef_12 = dist_12 * interval_12;
                    coef_21 = dist_21 * interval_21;
                    coef_22 = dist_22 * interval_22;

                    C_11 = coef_12 * coef_21 * coef_22;
                    C_12 = coef_11 * coef_21 * coef_22;
                    C_21 = coef_11 * coef_12 * coef_22;
                    C_22 = coef_11 * coef_12 * coef_21;

                    coef_sum = C_11 + C_12 + C_21 + C_22;
                    
                    ts_11 = obj.TimestampMap(j_floor + 1, i_floor + 1);
                    ts_12 = obj.TimestampMap(j_floor + 1, i_floor + 2);
                    ts_21 = obj.TimestampMap(j_floor + 2, i_floor + 1);
                    ts_22 = obj.TimestampMap(j_floor + 2, i_floor + 2);

                    thr = (ts_11*C_11 + ts_12 * C_12 + ts_21 * C_21 + ts_22 * C_22) / coef_sum;

                end
            end
        end

        function mainDisplay(obj)
            obj.displayFilter();
        end
        
    end
end

