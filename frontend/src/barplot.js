import { useEffect, useMemo, useRef, useState, useLayoutEffect } from "react";
import * as d3 from "d3";

const MARGIN = { top: 30, right: 30, bottom: 50, left: 50 };


const useDimensions = (targetRef) => {
    const getDimensions = () => {
        return {
            width: targetRef.current ? targetRef.current.offsetWidth : 0,
            height: targetRef.current ? targetRef.current.offsetHeight : 0
        };
    };

    const [dimensions, setDimensions] = useState(getDimensions);

    const handleResize = () => {
        setDimensions(getDimensions());
    };

    useEffect(() => {
        window.addEventListener("resize", handleResize);
        return () => window.removeEventListener("resize", handleResize);
    }, []);

    useLayoutEffect(() => {
        handleResize();
    }, []);

    return dimensions;
}


const getStackedBarplotData = (featureValues, selectedPoints, selectedGroupName) => {
    const uniqueValues = [...new Set(featureValues)];
    var selectedCounts = {};
    var otherCounts = {};
    uniqueValues.forEach((value) => {
        selectedCounts[value] = 0;
        otherCounts[value] = 0;
    });

    var isSelected = new Array(featureValues.length).fill(false);
    selectedPoints.forEach((point) => {
        isSelected[point] = true;
    });

    for (var i = 0; i < featureValues.length; i++) {
        var value = featureValues[i];
        if (isSelected[i]) {
            selectedCounts[value] = selectedCounts[value] + 1;
        } else {
            otherCounts[value] = otherCounts[value] + 1;
        }
    }

    // normalize by dividing by the number of selected/other points (selectedPoints.length)
    uniqueValues.forEach((value) => {
        selectedCounts[value] = selectedCounts[value] / (selectedPoints.length == 0 ? 1 : selectedPoints.length);
        otherCounts[value] = otherCounts[value] / (featureValues.length - selectedPoints.length);
    });

    selectedCounts["x"] = selectedGroupName;
    otherCounts["x"] = "other";

    return [selectedCounts, otherCounts];
}

const getDodgedBarplotData = (featureValues, selectedPoints, selectedGroupName) => {
    const uniqueValues = [...new Set(featureValues)];
    var selectedCounts = {};
    var otherCounts = {};
    uniqueValues.forEach((value) => {
        selectedCounts[value] = 0;
        otherCounts[value] = 0;
    });

    var isSelected = new Array(featureValues.length).fill(false);
    selectedPoints.forEach((point) => {
        isSelected[point] = true;
    });

    for (var i = 0; i < featureValues.length; i++) {
        var value = featureValues[i];
        if (isSelected[i]) {
            selectedCounts[value] = selectedCounts[value] + 1;
        } else {
            otherCounts[value] = otherCounts[value] + 1;
        }
    }

    // normalize by dividing by the number of selected/other points (selectedPoints.length)
    uniqueValues.forEach((value) => {
        selectedCounts[value] = selectedCounts[value] / (selectedPoints.length == 0 ? 1 : selectedPoints.length);
        otherCounts[value] = otherCounts[value] / (featureValues.length - selectedPoints.length);
    });

    const result = {};
    result[selectedGroupName] = selectedCounts;
    result["other"] = otherCounts;

    return result;
}


export const DodgedBarplot = ({
    featureValues,
    xlabel,
    selectedPoints,
    selectedGroupName,
    colorMap
}) => {
    // TODO only show the first 8-10 categories if there are more.
    // bounds = area inside the graph axis = calculated by substracting the margins
    const axesRef = useRef(null);
    const targetRef = useRef(null);

    const dimensions = useDimensions(targetRef);
    const width = dimensions.width;
    const height = dimensions.height;
    const boundsWidth = width - MARGIN.right - MARGIN.left;
    const boundsHeight = height - MARGIN.top - MARGIN.bottom;

    const data = getDodgedBarplotData(featureValues, selectedPoints, selectedGroupName);
    const categories = colorMap["ticks"];
    const groups = Object.keys(data);

    var max_y_value = Math.max(Math.max(...Object.values(data[selectedGroupName])), Math.max(...Object.values(data["other"])));
    max_y_value = Math.ceil(max_y_value * 10) / 10;

    // Y axis
    const yScale = useMemo(() => {
        return d3
            .scaleLinear()
            .domain([0, max_y_value])
            .range([boundsHeight, 0]);
    }, [data, height]);

    // X axis
    const xScale = useMemo(() => {
        return d3
            .scaleBand()
            .domain(categories)
            .range([0, boundsWidth])
            .padding(0.2);
    }, [data, width]);

    // Color Scale
    var colorScale = d3
        .scaleOrdinal()
        .domain(categories)
        .range(colorMap["colors"]);

    useEffect(() => {
        const svgElement = d3.select(axesRef.current);
        svgElement.selectAll("*").remove();
        const xAxisGenerator = d3.axisBottom(xScale);
        svgElement
            .append("g")
            .attr("transform", "translate(0," + boundsHeight + ")")
            .call(xAxisGenerator);

        svgElement.selectAll(".domain,.tick>line")
            .remove();

        const yAxisGenerator = d3.axisLeft(yScale);
        svgElement.append("g").call(yAxisGenerator);

        svgElement.append("text")
            .attr("y", boundsHeight + 40)
            .attr("x", boundsWidth / 2)
            .style("text-anchor", "middle")
            .text(xlabel);
    }, [xScale, yScale, boundsHeight]);

    const barwidth = 10 / 21 * xScale.bandwidth();
    const bar_spacing = 1 / 21 * xScale.bandwidth();

    const rectangles = groups.map((group, i) => {
        return (
            <g key={i}>
                {categories.map((category, j) => {
                    if (boundsHeight - yScale(data[group][category]) < 0) {
                        return null;
                    } else {
                        return (
                            <rect
                                key={j}
                                x={i == 0 ? xScale(category) : (xScale(category) + xScale.bandwidth() / 2) + bar_spacing}
                                y={yScale(data[group][category])}
                                height={boundsHeight - yScale(data[group][category])}
                                width={barwidth}
                                fill={colorScale(category)}
                                opacity={group == "other" ? 0.5 : 0.9}
                            ></rect>
                        );
                    }
                })}
            </g>
        );
    });

    return (
        <div className="flex flex-wrap items-center my-2 text-sm">
            <div className="flex flex-row w-full justify-end">
                {groups.map((group, i) => {
                    return <div key={i} className="flex flex-row items-center">
                        <svg width={12} height={12}><rect x={0} y={0} width={12}
                            height={12} fill={colorScale(categories[0])} opacity={group == "other" ? 0.5 : 0.9} /></svg>
                        <p className="text-sm ml-1 mr-4">{group}</p>
                    </div>;
                })
                }
            </div>
            <div className='flex flex-wrap items-center justify-left text-sm min-w-full max-w-full h-[240px]'
                ref={targetRef}>
                <svg width={width} height={height}>
                    <g
                        width={boundsWidth}
                        height={boundsHeight}
                        transform={`translate(${[MARGIN.left, MARGIN.top].join(",")})`}
                    >
                        {rectangles}
                    </g>
                    <g
                        width={boundsWidth}
                        height={boundsHeight}
                        ref={axesRef}
                        transform={`translate(${[MARGIN.left, MARGIN.top].join(",")})`}
                    />
                </svg>
            </div>
        </div>
    );

}


export const StackedBarplot = ({
    featureValues,
    xlabel,
    selectedPoints,
    selectedGroupName,
    colorMap
}) => {
    // bounds = area inside the graph axis = calculated by substracting the margins
    const axesRef = useRef(null);
    const targetRef = useRef(null);

    const dimensions = useDimensions(targetRef);
    const width = dimensions.width;
    const height = dimensions.height;
    const boundsWidth = width - MARGIN.right - MARGIN.left;
    const boundsHeight = height - MARGIN.top - MARGIN.bottom;


    const data = getStackedBarplotData(featureValues, selectedPoints, selectedGroupName);
    const allGroups = data.map((d) => String(d.x));

    const stackSeries = d3.stack().keys(colorMap["ticks"]);
    const series = stackSeries(data);

    // Y axis
    const yScale = useMemo(() => {
        return d3
            .scaleLinear()
            .domain([0, 1])
            .range([boundsHeight, 0]);
    }, [data, height]);

    // X axis
    const xScale = useMemo(() => {
        return d3
            .scaleBand()
            .domain(allGroups)
            .range([0, boundsWidth])
            .padding(0.05);
    }, [data, width]);

    // Color Scale
    var colorScale = d3
        .scaleOrdinal()
        .domain(colorMap["ticks"])
        .range(colorMap["colors"]);

    // Render the X and Y axis using d3.js, not react
    useEffect(() => {
        const svgElement = d3.select(axesRef.current);
        svgElement.selectAll("*").remove();
        const xAxisGenerator = d3.axisBottom(xScale);
        svgElement
            .append("g")
            .attr("transform", "translate(0," + boundsHeight + ")")
            .call(xAxisGenerator);

        svgElement.selectAll(".domain,.tick>line")
            .remove();

        const yAxisGenerator = d3.axisLeft(yScale);
        svgElement.append("g").call(yAxisGenerator);

        svgElement.append("text")
            .attr("y", boundsHeight + 40)
            .attr("x", boundsWidth / 2)
            .style("text-anchor", "middle")
            .text(xlabel);
    }, [xScale, yScale, boundsHeight]);

    const rectangles = series.map((subgroup, i) => {
        return (
            <g key={i}>
                {subgroup.map((group, j) => {
                    return (
                        <rect
                            key={j}
                            x={xScale(group.data.x)}
                            y={yScale(group[1])}
                            height={(yScale(group[0]) - yScale(group[1])) < 0 ? 0 : yScale(group[0]) - yScale(group[1])}
                            width={xScale.bandwidth()}
                            fill={colorScale(subgroup.key)}
                            opacity={0.9}
                        ></rect>
                    );
                })}
            </g>
        );
    });

    return (
        <div className='flex flex-wrap items-center justify-left text-sm min-w-full max-w-full h-[240px]'
            ref={targetRef}>
            <svg width={width} height={height}>
                <g
                    width={boundsWidth}
                    height={boundsHeight}
                    transform={`translate(${[MARGIN.left, MARGIN.top].join(",")})`}
                >
                    {rectangles}
                </g>
                <g
                    width={boundsWidth}
                    height={boundsHeight}
                    ref={axesRef}
                    transform={`translate(${[MARGIN.left, MARGIN.top].join(",")})`}
                />
            </svg>
        </div>
    );
};
