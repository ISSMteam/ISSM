function output_txt = plot_datatip(~, event_obj, x, y, field, datatip)

    pos = event_obj.Position;
    [~, idx] = min( (x - pos(1)).^2 + (y - pos(2)).^2 );
    
    if strcmp(datatip, 'vertex-field')
    output_txt = {
        ['X: ', num2str(pos(1))], ...
        ['Y: ', num2str(pos(2))], ...
        ['Vertex #: ', num2str(idx)], ...
        ['Field: ', num2str(field(idx))]
                 };
    else
        output_txt = {
            ['X: ', num2str(pos(1))], ...
            ['Y: ', num2str(pos(2))], ...
            ['Value: ', num2str(datatip(idx))]
                     };
    end

end
