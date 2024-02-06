import React, { useState } from 'react';
import { SettingsButton, ChevronRightButton, AsyncButton, DefaultButton } from "./buttons";


// Tab source: https://www.devwares.com/blog/how-to-create-react-tabs-with-tailwind-css/
const Tabs = ({ children }) => {
    const [activeTab, setActiveTab] = useState(children[0].props.label);
    const [visibility, setVisibility] = useState('visible');

    const handleClick = (e, newActiveTab) => {
        e.preventDefault();
        setActiveTab(newActiveTab);
    };

    const toggleVisibility = () => {
        if (visibility == "visible") setVisibility("hidden"); else setVisibility("visible");
    }

    if (visibility == "visible") {
        return (
            <div className="select-none right-0 w-1/4 min-w-[300px] h-screen max-h-screen bg-gray-100">
                <span className="relative">
                    <div className="absolute -left-[50px] top-2">
                        <ChevronRightButton onClick={toggleVisibility} />
                    </div>
                </span>

                <div className="flex border-b border-gray-300">
                    {children.map(child => (
                        <button
                            key={child.props.label}
                            className={`${activeTab === child.props.label ? 'border-b-2 border-purple-500 bg-gray-50' : ''
                                } flex-1 py-2 text-lg font-medium text-gray-900`}
                            onClick={e => handleClick(e, child.props.label)}
                        >
                            {child.props.label}
                        </button>
                    ))}
                </div>
                <div className="py-4">
                    {children.map(child => {
                        if (child.props.label === activeTab) {
                            return <div key={child.props.label}>{child.props.children}</div>;
                        }
                        return null;
                    })}
                </div>
            </div>
        );
    } else {
        return (
            <>
                <div className="fixed right-0 top-0 my-2">
                    <SettingsButton onClick={toggleVisibility} />
                </div>
            </>
        )
    }
};

const Tab = ({ label, children }) => {
    return (
        <div label={label} className="hidden">
            {children}
        </div>
    );
};
export { Tabs, Tab };