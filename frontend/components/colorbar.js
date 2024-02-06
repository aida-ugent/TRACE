export function Colorbar(props) {
    const { colormap, title, visibility } = props;
    let items, colors;

    if (colormap.constructor === 'Map') {
        items = Array.from(colormap.keys())
        colors = Array.from(colormap.values())
    } else {
        items = Object.keys(colormap);
        colors = Object.values(colormap);
    }

    // the labels will appear from largest to smallest
    if (parseFloat(items[0]) < parseFloat(items[1])) {
        items = items.reverse()
        colors = colors.reverse()
    }

    if (items.length > 11){
        // take every second element from items
        let newItems = items.filter((_, index) => index % 2 == 0)

        if (items.length % 2 == 0){
            // add the smallest element explicitly at the end
            newItems.push(items[items.length - 1])
        }
        items = newItems;
    }

    if (colors.length < 2) {
        console.log("Colorbar only got one color to create a linear gradient. Please specify at least two.")
        colors = colors.push("#a50026");
    }

    // the linear gradient is from top to bottom
    var bg_color = "linear-gradient(";
    for (let c of colors) {
        bg_color = bg_color + c + ", "
    }
    bg_color = bg_color.substring(0, bg_color.length - 2) + ")";

    return (
        <>
            {!(title == "none") ?
                (
                    <div className={visibility + " flex-col select-none h-fit w-fit max-w-[200px] overflow-auto bg-white/90 rounded-lg bg-clip-border mt-1 px-2 pt-1 pb-2"}>
                        <h4>{title}</h4>
                        <div className="flex flex-row min-h-[200px]">
                            <div className="flex select-none w-[20px] min-h-full rounded-sm bg-clip-border mt-2 mb-2"
                                style={{ background: bg_color }}>
                            </div>
                            <div className="flex min-w-fit flex-col justify-between">
                                {items.map((name) =>
                                    <p key={name}>&#x2012; {name}</p>
                                )}
                            </div>
                        </div >
                    </div>
                ) :
                (<></>)
            }
        </>
    )
};