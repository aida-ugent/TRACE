import React, { useState } from 'react';

// Tab source: https://www.devwares.com/blog/how-to-create-react-tabs-with-tailwind-css/
const Tabs = ({ children }) => {
    const [activeTab, setActiveTab] = useState(children[0].props.label);

    const handleClick = (e, newActiveTab) => {
        e.preventDefault();
        setActiveTab(newActiveTab);
    };


    return (
        // <div className="select-none w-full h-screen max-h-screen">
        <>
            <div className="flex flex-wrap content-normal border-b border-gray-300">
                {children.map(child => (
                    <button
                        key={child.props.label}
                        className={`${activeTab === child.props.label
                            ? 'border-b-2 border-purple-500 bg-gray-50'
                            : 'border-b-2 border-transparent'} 
                            min-w-fit grow py-2 pt-3 px-2 text-lg font-medium text-gray-900`}
                        onClick={e => handleClick(e, child.props.label)}
                    >
                        {child.props.label}
                    </button>
                ))}
            </div>
            <div className="py-4 px-5">
                {children.map(child => {
                    if (child.props.label === activeTab) {
                        return <div key={child.props.label}>{child.props.children}</div>;
                    } else {
                        return <div key={child.props.label} className='hidden'>{child.props.children}</div>
                    }
                })}
            </div>
        </>
    );
};

const Tab = ({ label, children }) => {
    return (
        <div label={label} className="hidden">
            {children}
        </div>
    );
};
export { Tabs, Tab };