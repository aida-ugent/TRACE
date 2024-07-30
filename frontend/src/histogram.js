import { useEffect, useMemo, useRef, useState, useLayoutEffect } from "react";
import * as d3 from "d3";

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


const MARGIN = { top: 10, right: 10, bottom: 75, left: 60 };
const BUCKET_NUMBER = 18;
const BUCKET_PADDING = 4;

// from https://www.react-graph-gallery.com/histogram
export const Histogram = ({ data, xlabel }) => {
    const targetRef = useRef(null);
    const dimensions = useDimensions(targetRef);
    const width = dimensions.width;
    const height = dimensions.height;

    const axesRef = useRef(null);
    const boundsWidth = width - MARGIN.right - MARGIN.left;
    const boundsHeight = height - MARGIN.top - MARGIN.bottom;

    const allGroupNames = data.map((group) => group.name);
    const allGroupColors = data.map((group) => group.color);
    const colorScale = d3
        .scaleOrdinal()
        .domain(allGroupNames)
        .range(allGroupColors);

    const xScale = useMemo(() => {
        const maxPerGroup = data.map((group) => Math.max(...group.values));
        const minPerGroup = data.map((group) => Math.min(...group.values));
        const max = Math.max(...maxPerGroup);
        const min = Math.min(...minPerGroup);
        return d3.scaleLinear().domain([min, max]).range([10, boundsWidth]).nice();
    }, [data, width]);

    const bucketGenerator = useMemo(() => {
        return d3
            .bin()
            .value((d) => d)
            .domain(xScale.domain())
            .thresholds(xScale.ticks(BUCKET_NUMBER)); // approximate number of bins (not exact!)
    }, [xScale]);

    // create categorical axis
    

    const groupBuckets = useMemo(() => {
        return data.map((group) => {
            return { name: group.name, buckets: bucketGenerator(group.values), size: group.values.length == 0 ? 1 : group.values.length };
        });
    }, [data]);

    const yScale = useMemo(() => {
        const max = Math.max(
            ...groupBuckets.map((group) => {
                return Math.max(...group.buckets.map((bucket) => (bucket?.length) / group.size));
            }
            )
        );
        return d3.scaleLinear().range([boundsHeight, 0]).domain([0, max]).nice();
    }, [data, height]);

    // Render the X axis using d3.js, not react
    useEffect(() => {
        const svgElement = d3.select(axesRef.current);
        svgElement.selectAll("*").remove();

        const xAxisGenerator = d3.axisBottom(xScale);

        svgElement
            .append("g")
            .attr("transform", "translate(0," + boundsHeight + ")")
            .call(xAxisGenerator);

        // rotate x-axis labels if they are too long
        const ticklabel_length = xAxisGenerator.scale().ticks().map(xScale.tickFormat())[0].length
        if (ticklabel_length > 3) {
            svgElement
                .selectAll("text")
                .style("text-anchor", "end")
                .attr("dx", "-.8em")
                .attr("dy", ".15em")
                .attr("transform", "rotate(-65)");
        }

        const yAxisGenerator = d3.axisLeft(yScale);
        svgElement.append("g").call(yAxisGenerator);


        // Add the text label for the Y axis
        svgElement.append("text")
            .attr("transform", "rotate(-90)")
            .attr("y", 0 - MARGIN.left)
            .attr("x", 0 - boundsHeight / 2)
            .attr("dy", "1em")
            .style("text-anchor", "middle")
            .text("relative frequency");

        svgElement.append("text")
            // .attr("transform", "translate(" + (width / 2) + " ," + (height) + ")")
            .attr("y", ticklabel_length > 3 ? boundsHeight + 70 : boundsHeight + 40)
            .attr("x", boundsWidth / 2)
            .style("text-anchor", "middle")
            .text(xlabel);

        svgElement.selectAll("text").attr("font-size", 14);
        svgElement.selectAll("path").style("stroke", "gray");
        svgElement.selectAll("line").style("stroke", "gray");
    }, [xScale, yScale, boundsHeight]);

    const allRects = groupBuckets.map((group, i) =>
        group.buckets.map((bucket, j) => {
            const { x0, x1 } = bucket;
            if (x0 == undefined || x1 == undefined) {
                return null;
            }
            //console.log(`group: ${group.name}, group size: ${group.size},  bucket: ${bucket.length}`);
            if (xScale(x1) - xScale(x0) - BUCKET_PADDING < 0) {
                console.log(`group: ${group.name}, x1: ${x1}, x0: ${x0}, xScale(x1): ${xScale(x1)}, xScale(x0): ${xScale(x0)}, BUCKET_PADDING: ${BUCKET_PADDING}`);
            }
            return (
                <rect
                    key={i + "_" + j}
                    fill={colorScale(group.name)}
                    opacity={0.6}
                    x={xScale(x0) + BUCKET_PADDING / 2}
                    width={xScale(x1) - xScale(x0) - BUCKET_PADDING}
                    y={yScale(bucket.length / group.size)}
                    height={boundsHeight - yScale(bucket.length / group.size)}
                />
            );
        })
    );

    return (
        <div className="flex flex-wrap items-center mb-2 justify-left text-sm">
            {groupBuckets.map(group => {
                return <>
                    <svg width={12} height={12}><rect x={0} y={0} width={12} height={12} fill={colorScale(group.name)} /></svg>
                    <p className="text-sm ml-1 mr-4">{group.name}</p>
                </>;
            })
            }
            <div className='flex flex-wrap items-center justify-left text-sm min-w-full max-w-full h-[220px]'
                ref={targetRef}>
                <svg width={width} height={height}>
                    <g
                        width={boundsWidth}
                        height={boundsHeight}
                        transform={`translate(${[MARGIN.left, MARGIN.top].join(",")})`}
                    >
                        {allRects}
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
};
