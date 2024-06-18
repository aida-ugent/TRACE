import React, { CSSProperties } from 'react';
import Select from 'react-select';

export default function GroupedSelect({ options, selected, onChange, menuPlacement = 'auto', isMulti = false, isDisabled = false, isClearable = false }) {

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


    if (isMulti) {
        //transform array of selected values into array with {'value': v, 'lavel': v} objects
        var selected_pairs = selected.map(v => {return { 'value': v, 'label': v }});
        return (
            <Select
                className='w-full text-slate-600 text-left'
                isClearable={isClearable}
                isSearchable={true}
                isDisabled={isDisabled}
                isMulti
                menuPlacement={menuPlacement}
                value={selected_pairs}
                onChange={(selection) => { onChange(selection); }}
                options={options}
                formatGroupLabel={formatGroupLabel}
                styles={{
                    control: (baseStyles, state) => ({
                        ...baseStyles,
                        whiteSpace: 'normal', //'pre-wrap',
                    }),
                }}
            />
        )
    } else {

        return (
            <Select
                className='w-full text-slate-600 text-left'
                isClearable={isClearable}
                isDisabled={isDisabled}
                isSearchable={true}
                menuPlacement={menuPlacement}
                value={{ 'value': selected, 'label': selected }}
                onChange={(selection) => { onChange(selection["value"]); }}
                options={options}
                formatGroupLabel={formatGroupLabel}
                styles={{
                    control: (baseStyles, state) => ({
                        ...baseStyles,
                        whiteSpace: 'normal', //'pre-wrap',
                    }),
                }}
            />
        )
    }
}