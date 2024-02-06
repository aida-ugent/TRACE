import { useState } from "react";

function LegendItem({ text, filter, unfilter, color }) {
    const [selected, setSelected] = useState(true);

    function handleClick() {
        if (selected) {
            filter(text)
            setSelected(false)
        } else {
            unfilter(text)
            setSelected(true)
        }
    }

    return (
        <li className="mt-1 flex cursor-pointer"
            key={text}
            onClick={handleClick}>
            <div className="flex-none w-4 h-4 mr-1 mt-1 rounded-full" style={{ 'backgroundColor': color }} />
            {selected ? (<p>{text}</p>) : (<p className="text-gray-400">{text}</p>)}
        </li>
    );
}

export default function Legend({ colormap, title, filter, unfilter, visibility }) {
    let items, colors;

    // print warning if there are less colors than items
    if (colormap.constructor === 'Map') {
        items = Array.from(colormap.keys())
        colors = Array.from(colormap.values())
    } else {
        items = Object.keys(colormap);
        colors = Object.values(colormap);
    }

    if (title === "none") {
        return (<></>)
    } else {
        return (
            <>
                <div className={visibility + " flex flex-col select-none w-fit max-w-[220px] overflow-auto bg-white/0 rounded-lg bg-clip-border mt-1 px-2 pt-1 pb-2"}>
                    <h4>{title}</h4>
                    <ul className=" min-w-fit h-fit">
                        {items.map((name, index) =>
                            <LegendItem key={name} text={name} filter={filter} unfilter={unfilter} color={colors[index]} />
                        )}
                    </ul>
                </div>
            </>
        );
    }
}
