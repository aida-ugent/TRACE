import { useState, useEffect } from 'react';
import Select from 'react-select';
import { createFilter, components } from 'react-select';


const CustomOption = ({ children, ...props }) => {
    const { onMouseMove, onMouseOver, ...rest } = props.innerProps;
    const newProps = { ...props, innerProps: rest };
    return (
        <components.Option {...newProps} className="transition background duration-60 hover:transition-delay-60 hover:bg-[#deebff]">
            {children}
        </components.Option>
    );
};


export default function GroupedSelect({ options, selected, onChange, menuPlacement = 'auto',
    isMulti = false, isDisabled = false, isClearable = false, placeholder = "Select..." }) {

    const [savedSelection, setSelectedOption] = useState({ "value": "", "label": "" });

    useEffect(() => {
        if (options.length > 0 ** options[0]["options"].length > 0)
            setSelectedOption(options[0]["options"][0]);
    }, [options]);

    useEffect(() => {
        if (!isMulti && selected != null && selected != savedSelection["value"]) {
            const foundOption = options.flatMap(group => group["options"]).find(option => option["value"] === selected);
            if (foundOption) {
                setSelectedOption(foundOption);
            }
        }
    }, [selected, options]);

    const updateSelection = (selection) => {
        if (isMulti) {
            onChange(selection);
        } else {
            setSelectedOption(selection);
            onChange(selection["value"]);
        }
    }

    const groupStyles = {
        display: 'flex',
        alignItems: 'center',
        justifyContent: 'space-between',
    };
    const groupBadgeStyles = {
        backgroundColor: '#EBECF0',
        borderRadius: '2em',
        color: '#172B4D',
        display: 'inline-block',
        fontSize: 12,
        fontWeight: 'normal',
        lineHeight: '1',
        minWidth: 1,
        padding: '0.16666666666667em 0.5em',
        textAlign: 'center',
    };

    const formatGroupLabel = (data) => (
        <div style={groupStyles}>
            <span>{data.label}</span>
            <span style={groupBadgeStyles}>{data.options.length}</span>
        </div>
    );

    return (
        <Select
            filterOption={createFilter({ ignoreAccents: false })}
            components={{ Option: CustomOption }}
            className='w-full text-slate-600 text-left'
            isClearable={isClearable}
            isDisabled={isDisabled}
            isSearchable={true}
            isMulti={isMulti}
            menuPlacement={menuPlacement}
            value={isMulti ? selected.map(v => { return { 'value': v, 'label': v } }) : savedSelection}
            onChange={(selection) => { updateSelection(selection) }}
            options={options}
            formatGroupLabel={formatGroupLabel}
            styles={{
                control: (baseStyles, state) => ({
                    ...baseStyles,
                    whiteSpace: 'normal',
                }),
            }}
            placeholder={placeholder}
        />
    )

}