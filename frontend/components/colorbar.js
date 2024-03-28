export function Colorbar(props) {
    const { colormap, title, visibility } = props;
    let items, colors;

    if (colormap instanceof Map) {
        items = Array.from(colormap.keys());
        colors = Array.from(colormap.values());
    } else {
        // copy ticks and color arrays
        items = Array.from(colormap.ticks);
        colors = Array.from(colormap.colors);
    }

    items = items.map(parseFloat);

    // Sort items in descending order
    if (items[0] < items[items.length - 1]) {
        items = items.reverse()
        colors = colors.reverse()
    }

    if (items.length > 11) {
        // Take every second element from items
        let newItems = items.filter((_, index) => index % 2 === 0);

        if (items.length % 2 === 0) {
            // Add the smallest element explicitly at the end
            newItems.push(items[items.length - 1]);
        }
        items = newItems;
    }

    // Calculate the smallest difference between two values
    let smallestDiff = Math.max(...items);
    for (let i = 0; i < items.length - 1; i++) {
        for (let j = i + 1; j < items.length; j++) {
            let diff = Math.abs(parseFloat(items[i]) - parseFloat(items[j]));
            if (diff < smallestDiff) {
                smallestDiff = diff;
            }
        }
    }

    // Format items with precision based on the smallest difference
    const match = smallestDiff.toString().match(/^[0\.]+/)
    const level = match ? match[0].length - 1 : 3;
    items = items.map((v) => v.toFixed(level));

    if (colors.length < 2) {
        console.log("Colorbar only got one color to create a linear gradient. Please specify at least two.");
        colors.push("#a50026");
    }

    // Create the linear gradient background color
    const bg_color = `linear-gradient(${colors.join(", ")})`;

    return (
        <>
            {title !== "none" && (
                <div className={`${visibility} flex-col select-none h-fit w-fit max-w-[200px] overflow-auto bg-white/90 rounded-lg bg-clip-border mt-1 px-2 pt-1 pb-2`}>
                    <h4>{title}</h4>
                    <div className="flex flex-row min-h-[200px]">
                        <div className="flex select-none w-[20px] min-h-full rounded-sm bg-clip-border mt-2 mb-2" style={{ background: bg_color }}></div>
                        <div className="flex min-w-fit flex-col justify-between">
                            {items.map((name, i) => (
                                <p key={`${i}_${name}`}>&#x2012; {name}</p>
                            ))}
                        </div>
                    </div>
                </div>
            )}
        </>
    );
}

